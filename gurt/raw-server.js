import net from "net";
import tls from "tls";
import fs from "fs";

export class GurtServer {
    constructor(options = {}) {
        this.host = options.host || "localhost";
        this.port = options.port || 4878;
        this.cert = options.cert;
        this.key = options.key;
        this.useTLS = !!(this.cert && this.key);
        this.server = null;
        this.connections = new Set();
        this.onRequest = null;
    }

    async listen() {
        return new Promise((resolve, reject) => {
            this.server = net.createServer((socket) => {
                this._handleConnection(socket);
            });

            this.server.listen(this.port, this.host, () => {
                resolve();
            });

            this.server.on("error", (err) => {
                reject(err);
            });
        });
    }

    _handleConnection(socket) {
        let buffer = "";
        let handshakeDone = false;
        let tlsSocket = null;
        let bodyLength = 0;
        const maxBodySize = 10 * 1024 * 1024;

        socket.setTimeout(20000);

        const handleData = (chunk) => {
            bodyLength += chunk.length;
            if (bodyLength > maxBodySize) {
                socket.destroy();
                return;
            }

            buffer += chunk.toString();

            if (!handshakeDone) {
                const handshakeEnd = buffer.indexOf("\r\n\r\n");
                if (handshakeEnd !== -1) {
                    const handshake = buffer.slice(0, handshakeEnd);
                    buffer = buffer.slice(handshakeEnd + 4);

                    if (this._handleHandshake(socket, handshake)) {
                        handshakeDone = true;

                        if (this.useTLS) {
                            try {
                                const secureContext = tls.createSecureContext({
                                    cert: fs.readFileSync(this.cert),
                                    key: fs.readFileSync(this.key)
                                });

                                socket.removeListener("data", handleData);

                                tlsSocket = new tls.TLSSocket(socket, {
                                    isServer: true,
                                    secureContext,
                                    ALPNProtocols: ["GURT/1.0"],
                                    requestCert: false,
                                    rejectUnauthorized: false
                                });

                                tlsSocket.once("secure", () => {
                                    this.connections.add(tlsSocket);
                                });

                                tlsSocket.on("data", handleData);
                                tlsSocket.on("error", (err) => {
                                    this.connections.delete(tlsSocket);
                                });
                                tlsSocket.on("close", () => {
                                    this.connections.delete(tlsSocket);
                                });
                            } catch (e) {
                                socket.destroy(e);
                                return;
                            }
                        } else {
                            this.connections.add(socket);
                        }
                    } else {
                        socket.end();
                    }
                }
            } else {
                this._handleRequest(tlsSocket || socket, buffer);
                buffer = "";
            }
        };

        socket.on("data", handleData);
        socket.on("error", (err) => {
            this.connections.delete(socket);
            if (tlsSocket) {
                this.connections.delete(tlsSocket);
            }
        });
        socket.on("close", () => {
            this.connections.delete(socket);
            if (tlsSocket) {
                this.connections.delete(tlsSocket);
            }
        });
        socket.on("timeout", () => {
            socket.destroy();
        });
    }

    _handleHandshake(socket, handshake) {
        const lines = handshake.split("\r\n");

        if (lines[0] === "HANDSHAKE / GURT/1.0.0") {

            const response = [
                "GURT/1.0.0 101 SWITCHING_PROTOCOLS",
                "gurt-version: 1.0.0",
                `encryption: ${this.useTLS ? "TLS/1.3" : "none"}`,
                "alpn: GURT/1.0",
                "server: GURT/1.0.0",
                `date: ${new Date().toUTCString()}`,
                "",
                ""
            ].join("\r\n");

            socket.write(response);
            return true;
        }

        socket.end();
        return false;
    }

    _handleRequest(socket, data) {
        const headerEnd = data.indexOf("\r\n\r\n");
        if (headerEnd === -1) {
            return;
        }

        const head = data.slice(0, headerEnd);
        const [requestLine, ...headerLines] = head.split("\r\n");
        const [method, path, version] = requestLine.split(" ");

        const headers = {};
        for (const h of headerLines) {
            const [k, v] = h.split(":");
            headers[k.trim().toLowerCase()] = v?.trim() ?? "";
        }

        const body = data.slice(headerEnd + 4);
        const contentLength = parseInt(headers["content-length"] || "0", 10);

        if (body.length < contentLength) {
            return;
        }

        const request = {
            method,
            path,
            headers,
            body: body.slice(0, contentLength)
        };

        if (this.onRequest) {
            this.onRequest(request, (response) => {
                this._sendResponse(socket, response);
            });
        } else {
        }
    }

    _sendResponse(socket, { code = 200, message = "OK", headers = {}, body = "" }) {
        const lines = [];
        lines.push(`GURT/1.0.0 ${code} ${message}`);
        headers = {
            "content-length": Buffer.byteLength(body),
            "server": "GURT/1.0.0",
            "date": new Date().toUTCString(),
            ...headers
        };
        for (const [k, v] of Object.entries(headers)) {
            lines.push(`${k}: ${v}`);
        }
        lines.push("", body);
        const responseData = lines.join("\r\n");
        socket.write(responseData);
    }

    close() {
        this.connections.forEach(conn => {
            conn.end();
        });
        this.connections.clear();
        if (this.server) {
            this.server.close();
        }
    }
}
