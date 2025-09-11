import { GurtApp } from "./gurt/app.js"

const app = new GurtApp({
    host: "mineshit.based",
    port: 4880,
    cert: "./mineshit.based.crt",
    key: "./mineshit.based.key"
})

app.use((req, res, next) => {
    console.log(`${req.method} ${req.path}`)
    next()
})

app.static("./root")


app.listen(4880, "0.0.0.0").then(() => {
    console.log("Server started on port 4880")
}).catch(console.error)

process.on("SIGINT", () => {
    console.log("Shutting down server...")
    app.close()
    process.exit(0)
})
