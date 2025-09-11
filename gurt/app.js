import { GurtServer } from "./raw-server.js";
import fs from "fs";
import path from "path";

export class GurtRequest {
    constructor(method, path, headers, body) {
        this.method = method;
        this.path = path;
        this.headers = headers;
        this.body = body;
        this.params = {};
        this.query = {};
        this._parseQuery();
    }

    _parseQuery() {
        const [pathname, queryString] = this.path.split("?");
        this.path = pathname;
        if (queryString) {
            for (const pair of queryString.split("&")) {
                const [key, value] = pair.split("=");
                this.query[decodeURIComponent(key)] = decodeURIComponent(value || "");
            }
        }
    }
}

export class GurtResponse {
    constructor(sendResponse) {
        this.sendResponse = sendResponse;
        this.code = 200;
        this.message = "OK";
        this.headers = {};
        this.body = "";
        this.sent = false;
    }

    status(code, message = "") {
        this.code = code;
        this.message = message;
        return this;
    }

    set(key, value) {
        this.headers[key] = value;
        return this;
    }

    send(body) {
        if (this.sent) return;
        this.body = body;
        this.sent = true;
        this.sendResponse(this);
        return this;
    }
    sendFile(filePath) {
        if (this.sent) return;
        try {
            const data = fs.readFileSync(filePath);
            if (!this.headers["content-type"]) {
                this.headers["content-type"] = this._getMimeType(filePath);
            }
            this.body = data;
            this.sent = true;
            this.sendResponse(this);
        } catch (err) {
            this.status(500, "Internal Server Error").json({ error: "Failed to read file" });
        }
        return this;
    }
    _getMimeType(filePath) {
        const ext = filePath.split(".").pop().toLowerCase();
        const mimeTypes = {
            html: "text/html",
            js: "application/javascript",
            css: "text/css",
            json: "application/json",
            png: "image/png",
            jpg: "image/jpeg",
            jpeg: "image/jpeg",
            gif: "image/gif",
            txt: "text/plain"
        };
        return mimeTypes[ext] || "application/octet-stream";
    }

    json(obj) {
        if (this.sent) return;
        this.headers["content-type"] = "application/json";
        this.body = JSON.stringify(obj);
        this.sent = true;
        this.sendResponse(this);
        return this;
    }
}

export class GurtApp {
    constructor(options = {}) {
        this.server = new GurtServer(options);
        this.middlewares = [];
        this.routes = {};
        this.server.onRequest = (req, sendResponse) => {
            this._handleRequest(req, sendResponse);
        };
    }

    use(middleware) {
        this.middlewares.push(middleware);
    }

    get(path, handler) {
        this._addRoute("GET", path, handler);
    }

    post(path, handler) {
        this._addRoute("POST", path, handler);
    }

    put(path, handler) {
        this._addRoute("PUT", path, handler);
    }

    delete(path, handler) {
        this._addRoute("DELETE", path, handler);
    }

    static(staticPath, options = {}) {
        const middleware = (req, res, next) => {
            if (req.method !== "GET") {
                return next();
            }
            
            let filePath = path.join(staticPath, req.path);
            
            if (req.path.endsWith("/")) {
                filePath = path.join(staticPath, req.path, options.index || "index.html");
            }
            
            try {
                const resolvedPath = path.resolve(filePath);
                const resolvedStatic = path.resolve(staticPath);
                
                if (!resolvedPath.startsWith(resolvedStatic)) {
                    return next();
                }
                
                if (fs.existsSync(resolvedPath) && fs.statSync(resolvedPath).isFile()) {
                    res.sendFile(resolvedPath);
                    return;
                }
            } catch (err) {
                
            }
            
            next();
        };
        
        this.use(middleware);
    }

    _addRoute(method, path, handler) {
        if (!this.routes[method]) {
            this.routes[method] = [];
        }
        this.routes[method].push({ path: this._pathToRegex(path), handler });
    }

    _pathToRegex(path) {
        const regex = path.replace(/:(\w+)/g, "(?<$1>[^/]+)");
        return new RegExp(`^${regex}$`);
    }

    async _handleRequest(rawReq, sendResponse) {
        const req = new GurtRequest(rawReq.method, rawReq.path, rawReq.headers, rawReq.body);
        const res = new GurtResponse(sendResponse);

        let index = 0;

        const next = async () => {
            if (index < this.middlewares.length) {
                await this.middlewares[index++](req, res, next);
            } else {
                await this._handleRoute(req, res);
            }
        };

        try {
            await next();
        } catch (err) {
            if (!res.sent) {
                res.status(500, "Internal Server Error").json({ error: err.message });
            }
        }

        if (!res.sent) {
            res.json({ message: "OK" });
        }
    }

    async _handleRoute(req, res) {
        const routes = this.routes[req.method] || [];
        for (const route of routes) {
            const match = req.path.match(route.path);
            if (match) {
                req.params = match.groups || {};
                await route.handler(req, res);
                return;
            }
        }
        res.status(404, "Not Found").json({ error: "Not found" });
    }

    listen(port, host = "localhost") {
        this.server.port = port;
        this.server.host = host;
        return this.server.listen();
    }

    close() {
        this.server.close();
    }
}
