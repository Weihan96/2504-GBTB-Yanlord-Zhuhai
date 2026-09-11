// Adapter to the installed public bridge's framed protocol; no server changes.
import net from 'node:net';
const [port, command, params = '{}'] = process.argv.slice(2);
if (!port || !command) throw new Error('Usage: bun call_review_bonsai_bridge.mjs PORT COMMAND PARAMS_JSON');
const payload = Buffer.from(JSON.stringify({id: 'review-storage-pilot', command, params: JSON.parse(params)}));
const header = Buffer.alloc(4); header.writeUInt32BE(payload.length);
await new Promise((resolve, reject) => {
  const socket = net.createConnection({host: '127.0.0.1', port: Number(port)});
  let data = Buffer.alloc(0);
  socket.setTimeout(600000);
  socket.on('connect', () => socket.write(Buffer.concat([header, payload])));
  socket.on('error', reject);
  socket.on('timeout', () => {socket.destroy(); reject(new Error('Bridge timeout; inspect result before retry'));});
  socket.on('data', chunk => {
    data = Buffer.concat([data, chunk]);
    if (data.length < 4 || data.length < 4 + data.readUInt32BE(0)) return;
    const result = JSON.parse(data.subarray(4, 4 + data.readUInt32BE(0)).toString());
    socket.end(); console.log(JSON.stringify(result));
    if (!result.success || result.result?.success === false) process.exitCode = 1;
    resolve();
  });
});
