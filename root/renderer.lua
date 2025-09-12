--!optimize 2
--!native
local canvas = gurt.select("#poly-test")
local ctx = canvas:withContext("2d")
local Renderer = { 
    meshes = {}, 
    camera = { pos = { 0, -30, -3 }, rot = { 0, 0, 0 }, fov = 200, cx = 400, cy = 300 },
    sun = { dir = { 0.5, -1, 0.3 }, ambient = 0.3, diffuse = 0.7 },
    hitFace = nil,
    worldChunks = {}
}

local keys = {}

local function rotXYZ(x, y, z, rx, ry, rz)
    local cy = math.cos(ry)
    local sy = math.sin(ry)
    local cx = math.cos(rx)
    local sx = math.sin(rx)
    local cz = math.cos(rz)
    local sz = math.sin(rz)
    local y1 = y * cx - z * sx
    local z1 = y * sx + z * cx
    local x2 = x * cy + z1 * sy
    local z2 = -x * sy + z1 * cy
    local x3 = x2 * cz - y1 * sz
    local y3 = x2 * sz + y1 * cz
    return x3, y3, z2
end

local function toView(cam, x, y, z)
    local vx = x - cam.pos[1]
    local vy = -y - cam.pos[2]
    local vz = z - cam.pos[3]
    local cx = math.cos(cam.rot[1])
    local sx = math.sin(cam.rot[1])
    local cy = math.cos(cam.rot[2])
    local sy = math.sin(cam.rot[2])
    local x1 = vx * cy + vz * sy
    local z1 = -vx * sy + vz * cy
    local y2 = vy * cx - z1 * sx
    local z2 = vy * sx + z1 * cx
    return x1, y2, z2
end

local function projectView(cam, vx, vy, vz)
    local s = cam.fov / vz
    return cam.cx + vx * s, cam.cy + vy * s, vz
end

local function normalize(x, y, z)
    local len = math.sqrt(x*x + y*y + z*z)
    if len == 0 then return 0, 0, 0 end
    return x/len, y/len, z/len
end

local function cross(ax, ay, az, bx, by, bz)
    return ay*bz - az*by, az*bx - ax*bz, ax*by - ay*bx
end

local function dot(ax, ay, az, bx, by, bz)
    return ax*bx + ay*by + az*bz
end

local function getRayFromCamera(cam)
    local cy = math.cos(cam.rot[2])
    local sy = math.sin(cam.rot[2])
    local cx = math.cos(cam.rot[1])
    local sx = math.sin(cam.rot[1])
    local dx = -sy * cx
    local dy = -sx
    local dz =  cy * cx
    local inv = 1 / math.sqrt(dx*dx + dy*dy + dz*dz)
    return { dx*inv, dy*inv, dz*inv }
end



local function calculateLighting(nx, ny, nz, sun)
    local sunx, suny, sunz = normalize(sun.dir[1], sun.dir[2], sun.dir[3])
    local dotProduct = -dot(nx, ny, nz, sunx, suny, sunz)
    local intensity = math.max(0, dotProduct) * sun.diffuse + sun.ambient
    return math.min(1, intensity)
end

local function rgbToHex(r, g, b)
    local rInt = math.floor(r * 255)
    local gInt = math.floor(g * 255)
    local bInt = math.floor(b * 255)
    return string.format("#%02x%02x%02x", rInt, gInt, bInt)
end

local function applyLighting(baseColor, intensity)
    if not baseColor or baseColor:sub(1,1) ~= "#" then
        return rgbToHex(intensity, intensity, intensity)
    end
    
    local rHex = tonumber(baseColor:sub(2,3), 16)
    local gHex = tonumber(baseColor:sub(4,5), 16)
    local bHex = tonumber(baseColor:sub(6,7), 16)
    
    if not rHex or not gHex or not bHex then
        return rgbToHex(intensity, intensity, intensity)
    end
    
    local r = rHex / 255
    local g = gHex / 255
    local b = bHex / 255
    
    r = r * intensity
    g = g * intensity
    b = b * intensity
    
    return rgbToHex(r, g, b)
end

local function lerp(a, b, t)
    return a + t * (b - a)
end

local function fade(t)
    return t * t * t * (t * (t * 6 - 15) + 10)
end

local function grad(hash, x, y, z)
    local h = hash % 16
    local u = h < 8 and x or y
    local v = h < 4 and y or (h == 12 or h == 14) and x or z
    return ((h % 2) == 0 and u or -u) + ((h % 4) < 2 and v or -v)
end

local perm = {}
for i = 0, 255 do
    perm[i] = math.random(0, 255)
    perm[i + 256] = perm[i]
end

local function noise3d(x, y, z)
    local X = math.floor(x) % 256
    local Y = math.floor(y) % 256
    local Z = math.floor(z) % 256
    
    x = x - math.floor(x)
    y = y - math.floor(y)
    z = z - math.floor(z)
    
    local u = fade(x)
    local v = fade(y)
    local w = fade(z)
    
    local A = perm[X] + Y
    local AA = perm[A] + Z
    local AB = perm[A + 1] + Z
    local B = perm[X + 1] + Y
    local BA = perm[B] + Z
    local BB = perm[B + 1] + Z
    
    return lerp(w, lerp(v, lerp(u, grad(perm[AA], x, y, z),
                                   grad(perm[BA], x - 1, y, z)),
                          lerp(u, grad(perm[AB], x, y - 1, z),
                                   grad(perm[BB], x - 1, y - 1, z))),
                   lerp(v, lerp(u, grad(perm[AA + 1], x, y, z - 1),
                                   grad(perm[BA + 1], x - 1, y, z - 1)),
                          lerp(u, grad(perm[AB + 1], x, y - 1, z - 1),
                                   grad(perm[BB + 1], x - 1, y - 1, z - 1))))
end

if not math.noise then
    math.noise = function(x, y, z)
        return noise3d(x or 0, y or 0, z or 0)
    end
end
local CHUNK = 8
local SCREEN_W, SCREEN_H = 800, 600
local MATERIAL_COLOR = {
    [1] = "#00aa00",
    [2] = "#777777"
}
local function clipPlane3D(poly, a,b,c,d)
    local out = {}
    local p0 = poly[#poly]
    local s0 = a*p0[1] + b*p0[2] + c*p0[3] + d
    local in0 = s0 >= 0
    for i=1,#poly do
        local p1 = poly[i]
        local s1 = a*p1[1] + b*p1[2] + c*p1[3] + d
        local in1 = s1 >= 0
        if in1 ~= in0 then
            local t = s0 / (s0 - s1)
            out[#out+1] = { p0[1] + (p1[1]-p0[1])*t, p0[2] + (p1[2]-p0[2])*t, p0[3] + (p1[3]-p0[3])*t }
        end
        if in1 then out[#out+1] = p1 end
        p0, s0, in0 = p1, s1, in1
    end
    return out
end

local function clipFrustumView(poly, cam)
    return clipPlane3D(poly, 0, 0, 1, -0.01)
end

local function drawPoly(points, color)
    if #points < 3 then return end
    ctx:beginPath()
    ctx:moveTo(points[1][1], points[1][2])
    for i = 2, #points do
        ctx:lineTo(points[i][1], points[i][2])
    end
    ctx:closePath()
    ctx:setFillStyle(color)
    ctx:fill()
end

local function strokePoly(points, color, lineWidth)
    if #points < 2 then return end
    ctx:beginPath()
    ctx:moveTo(points[1][1], points[1][2])
    for i = 2, #points do
        ctx:lineTo(points[i][1], points[i][2])
    end
    ctx:closePath()
    if color then ctx:setStrokeStyle(color) end
    if lineWidth then ctx:setLineWidth(lineWidth) end
    ctx:stroke()
end

function Renderer:addMesh(verts, faces, opts)
    local m = {
        vertices = verts,
        faces = faces,
        position = opts and opts.position or { 0, 0, 0 },
        rotation = opts and opts.rotation or { 0, 0, 0 },
        scale = opts and opts.scale or { 1, 1, 1 },
        color = opts and opts.color or nil
    }
    table.insert(self.meshes, m)
    return m
end

function Renderer:clear()
    self.meshes = {}
end

function Renderer:raycast()
    local o = table.clone(self.camera.pos)
    o[2] = -o[2]
    local d = getRayFromCamera(self.camera)
    local x, y, z = o[1], o[2], o[3]
    local dx, dy, dz = d[1], d[2], d[3]

    local stepX = dx > 0 and 1 or dx < 0 and -1 or 0
    local stepY = dy > 0 and 1 or dy < 0 and -1 or 0
    local stepZ = dz > 0 and 1 or dz < 0 and -1 or 0

    local tDeltaX = stepX ~= 0 and math.abs(1 / dx) or math.huge
    local tDeltaY = stepY ~= 0 and math.abs(1 / dy) or math.huge
    local tDeltaZ = stepZ ~= 0 and math.abs(1 / dz) or math.huge

    local mapX = math.floor(x)
    local mapY = math.floor(y)
    local mapZ = math.floor(z)

    local function tmax(pos, cell, step, delta)
        if step > 0 then return (cell + 1 - pos) * delta
        elseif step < 0 then return (pos - cell) * delta
        else return math.huge end
    end

    local tMaxX = tmax(x, mapX, stepX, tDeltaX)
    local tMaxY = tmax(y, mapY, stepY, tDeltaY)
    local tMaxZ = tmax(z, mapZ, stepZ, tDeltaZ)

    local nx, ny, nz = 0, 0, 0
    for _ = 1, 128 do
        if self:isBlockSolid(mapX, mapY, mapZ) then
            self.hitFace = { x = mapX, y = mapY, z = mapZ, nx = nx, ny = ny, nz = nz }
            return
        end
        if tMaxX < tMaxY and tMaxX < tMaxZ then
            mapX = mapX + stepX
            tMaxX = tMaxX + tDeltaX
            nx, ny, nz = -stepX, 0, 0
        elseif tMaxY < tMaxZ then
            mapY = mapY + stepY
            tMaxY = tMaxY + tDeltaY
            nx, ny, nz = 0, -stepY, 0
        else
            mapZ = mapZ + stepZ
            tMaxZ = tMaxZ + tDeltaZ
            nx, ny, nz = 0, 0, -stepZ
        end
    end
    self.hitFace = nil
end


local WORLD_Y = 32

function Renderer:isBlockSolid(wx, wy, wz)
    if wy < 0 or wy >= WORLD_Y then return false end
    local cx = math.floor(wx / CHUNK)
    local cz = math.floor(wz / CHUNK)
    local chunk = self.worldChunks and self.worldChunks[cx.."|"..cz]
    if not chunk then return false end
    local lx = wx - cx*CHUNK
    local lz = wz - cz*CHUNK
    local idx = wy*CHUNK*CHUNK + lz*CHUNK + lx + 1
    return chunk.data[idx] ~= 0
end


function Renderer:render()
    ctx:clearRect(0, 0, SCREEN_W, SCREEN_H)
    ctx:fillRect(0, 0, SCREEN_W, SCREEN_H, "#87CEEB")
    local all = {}
    for _, m in ipairs(self.meshes) do
        local vv, wv = {}, {}
        for i, v in ipairs(m.vertices) do
            local x = v[1] * m.scale[1]
            local y = v[2] * m.scale[2]
            local z = v[3] * m.scale[3]
            local x1, y1, z1 = rotXYZ(x, y, z, m.rotation[1], m.rotation[2], m.rotation[3])
            local wx = x1 + m.position[1]
            local wy = y1 + m.position[2]
            local wz = z1 + m.position[3]
            wv[i] = { wx, wy, wz }
            local vx, vy, vz = toView(self.camera, wx, wy, wz)
            vv[i] = { vx, vy, vz }
        end
        for _, f in ipairs(m.faces) do
            local i1, i2, i3, i4 = f[1], f[2], f[3], f[4]
            local wa, wb, wc = wv[i1], wv[i2], wv[i3]
            local v1x, v1y, v1z = wb[1] - wa[1], wb[2] - wa[2], wb[3] - wa[3]
            local v2x, v2y, v2z = wc[1] - wa[1], wc[2] - wa[2], wc[3] - wa[3]
            local nx, ny, nz = cross(v1x, v1y, v1z, v2x, v2y, v2z)
            nx, ny, nz = normalize(nx, ny, nz)

            local cxw, cyw, czw
            if i4 then
                local wd = wv[i4]
                cxw = 0.25*(wa[1]+wb[1]+wc[1]+wd[1])
                cyw = 0.25*(wa[2]+wb[2]+wc[2]+wd[2])
                czw = 0.25*(wa[3]+wb[3]+wc[3]+wd[3])
            else
                cxw = (wa[1]+wb[1]+wc[1])*(1/3)
                cyw = (wa[2]+wb[2]+wc[2])*(1/3)
                czw = (wa[3]+wb[3]+wc[3])*(1/3)
            end
            local camwx = self.camera.pos[1]
            local camwy = -self.camera.pos[2]
            local camwz = self.camera.pos[3]
            local ddx = camwx - cxw
            local ddy = camwy - cyw
            local ddz = camwz - czw
            if nx*ddx + ny*ddy + nz*ddz > 0 then
                local intensity = calculateLighting(nx, ny, nz, self.sun)
                local baseColor = m.color or "#888888"
                local faceColor = (f.color or f[5])
                local baseColor = faceColor or m.color or "#888888"
                local litColor = applyLighting(baseColor, intensity)

                local poly = {}
                poly[#poly+1] = { vv[i1][1], vv[i1][2], vv[i1][3] }
                poly[#poly+1] = { vv[i2][1], vv[i2][2], vv[i2][3] }
                poly[#poly+1] = { vv[i3][1], vv[i3][2], vv[i3][3] }
                if i4 then poly[#poly+1] = { vv[i4][1], vv[i4][2], vv[i4][3] } end
                poly = clipFrustumView(poly, self.camera)
                if #poly >= 3 then
                    local sp = {}
                    local zsum = 0
                    for i=1,#poly do
                        local px, py, pz = projectView(self.camera, poly[i][1], poly[i][2], poly[i][3])
                        sp[i] = { px, py, pz }
                        zsum = zsum + pz
                    end
                    local key = zsum / #sp
                    all[#all+1] = { sp, litColor, key }
                end
            end
        end
    end
    table.sort(all, function(a,b) return a[3] > b[3] end)
    for i=1,#all do
        local t = all[i]
        drawPoly(t[1], t[2])
    end
    
    if self.hitFace then
    local h = self.hitFace
    local q
    local hnx = h.nx or 0
    local hny = h.ny or 0
    local hnz = h.nz or 0
    if hnx ~= 0 then
        local X = hnx > 0 and h.x + 1 or h.x
        q = { {X,h.y,h.z}, {X,h.y+1,h.z}, {X,h.y+1,h.z+1}, {X,h.y,h.z+1} }
    elseif hny ~= 0 then
        local Y = hny > 0 and h.y + 1 or h.y
        q = { {h.x,Y,h.z}, {h.x+1,Y,h.z}, {h.x+1,Y,h.z+1}, {h.x,Y,h.z+1} }
    else
        local Z = hnz > 0 and h.z + 1 or h.z
        q = { {h.x,h.y,Z}, {h.x+1,h.y,Z}, {h.x+1,h.y+1,Z}, {h.x,h.y+1,Z} }
    end

    local poly = {}
    for i=1,4 do
        local vx,vy,vz = toView(self.camera, q[i][1], q[i][2], q[i][3])
        poly[i] = {vx,vy,vz}
    end
    poly = clipFrustumView(poly, self.camera)
    if #poly >= 3 then
        local sp = {}
        for i=1,#poly do
            local px,py,pz = projectView(self.camera, poly[i][1], poly[i][2], poly[i][3])
            sp[i] = {px,py,pz}
        end
    strokePoly(sp, "#ffff00", 2)
    end
end

end
local loadedMeshes = {}


local cubeVerts = {
    { -0.5, -0.5, -0.5 },
    {  0.5, -0.5, -0.5 },
    {  0.5,  0.5, -0.5 },
    { -0.5,  0.5, -0.5 },
    { -0.5, -0.5,  0.5 },
    {  0.5, -0.5,  0.5 },
    {  0.5,  0.5,  0.5 },
    { -0.5,  0.5,  0.5 }
}

local cubeFaces = {
    {1,4,3,2},
    {5,6,7,8},
    {1,5,8,4},
    {2,3,7,6},
    {4,8,7,3},
    {1,2,6,5}
}


local cube = Renderer:addMesh(cubeVerts, cubeFaces, { position = { 0, 2, 0 }, rotation = { 0, 0, 0 }, scale = { 1, 1, 1 }, color = "#ff0000" })

local SCALE = 1
local currentBlockId = 1

local World = {chunks = {}}
local function k(cx,cz) return cx.."|"..cz end
local function idx(x,y,z) return y*CHUNK*CHUNK + z*CHUNK + x + 1 end

function World:getChunk(cx,cz) return self.chunks[k(cx,cz)] end
function World:newChunk(cx,cz)
    local c = {cx=cx, cz=cz, data={}, verts=nil, faces=nil, dirty=false}
    local n = CHUNK*WORLD_Y*CHUNK
    for i=1,n do c.data[i]=0 end
    self.chunks[k(cx,cz)] = c
    Renderer.worldChunks[k(cx,cz)] = c
    return c
end
function World:placeBlock(wx, wy, wz, val)
    if wy < 0 or wy >= WORLD_Y then return false end
    local cx = math.floor(wx / CHUNK)
    local cz = math.floor(wz / CHUNK)
    local chunk = self:getChunk(cx, cz)
    if not chunk then return false end
    local lx = wx - cx * CHUNK
    local lz = wz - cz * CHUNK
    chunk.data[idx(lx, wy, lz)] = val or 1
    chunk.dirty = true
    self:remeshChunk(cx, cz)
    if lx == 0 then self:remeshChunk(cx-1, cz) end
    if lx == CHUNK-1 then self:remeshChunk(cx+1, cz) end
    if lz == 0 then self:remeshChunk(cx, cz-1) end
    if lz == CHUNK-1 then self:remeshChunk(cx, cz+1) end
    self:saveChunkToCrumb(cx, cz)
    return true
end

function World:removeBlock(wx, wy, wz)
    if wy < 0 or wy >= WORLD_Y then return false end
    local cx = math.floor(wx / CHUNK)
    local cz = math.floor(wz / CHUNK)
    local chunk = self:getChunk(cx, cz)
    if not chunk then return false end
    local lx = wx - cx * CHUNK
    local lz = wz - cz * CHUNK
    chunk.data[idx(lx, wy, lz)] = 0
    chunk.dirty = true
    self:remeshChunk(cx, cz)
    if lx == 0 then self:remeshChunk(cx-1, cz) end
    if lx == CHUNK-1 then self:remeshChunk(cx+1, cz) end
    if lz == 0 then self:remeshChunk(cx, cz-1) end
    if lz == CHUNK-1 then self:remeshChunk(cx, cz+1) end
    self:saveChunkToCrumb(cx, cz)
    return true
end

function World:remeshChunk(cx, cz)
    local verts, faces = self:meshChunk(cx, cz)
    local key = k(cx, cz)
    local mesh = loadedMeshes[key]
    if mesh then
        mesh.vertices = verts
        mesh.faces = faces
    end
end

function World:setBlock(wx,wy,wz,val)
    if wy<0 or wy>=WORLD_Y then return end
    local cx = math.floor(wx/CHUNK)
    local cz = math.floor(wz/CHUNK)
    local lx = wx - cx*CHUNK
    local lz = wz - cz*CHUNK
    local c = self:getChunk(cx,cz) or self:newChunk(cx,cz)
    c.data[idx(lx,wy,lz)] = val
end
function World:isSolid(wx,wy,wz)
    if wy<0 then return true end
    if wy>=WORLD_Y then return false end
    local cx = math.floor(wx/CHUNK)
    local cz = math.floor(wz/CHUNK)
    local c = self:getChunk(cx,cz)
    if not c then return false end
    local lx = wx - cx*CHUNK
    local lz = wz - cz*CHUNK
    return c.data[idx(lx,wy,lz)] ~= 0
end

function World:isSolidForMeshing(wx,wy,wz)
    if wy < 0 then return true end
    if wy >= WORLD_Y then return false end
    local cx = math.floor(wx/CHUNK)
    local cz = math.floor(wz/CHUNK)
    local c = self:getChunk(cx,cz)
    if not c then return false end
    local lx = wx - cx*CHUNK
    local lz = wz - cz*CHUNK
    return c.data[idx(lx,wy,lz)] ~= 0
end

function World:getBlock(wx,wy,wz)
    if wy < 0 or wy >= WORLD_Y then return 0 end
    local cx = math.floor(wx/CHUNK)
    local cz = math.floor(wz/CHUNK)
    local c = self:getChunk(cx,cz)
    if not c then return 0 end
    local lx = wx - cx*CHUNK
    local lz = wz - cz*CHUNK
    return c.data[idx(lx,wy,lz)] or 0
end

function World:saveToCrumb()
    local nPer = CHUNK*WORLD_Y*CHUNK
    local list = {}
    for _, c in pairs(self.chunks) do list[#list+1] = c end
    local function pack32(n)
        local u = n % 4294967296
        local b1 = u % 256; u = math.floor(u/256)
        local b2 = u % 256; u = math.floor(u/256)
        local b3 = u % 256; u = math.floor(u/256)
        local b4 = u % 256
        return string.char(b1,b2,b3,b4)
    end
    local parts = {}
    parts[#parts+1] = "WLD1"..pack32(CHUNK)..pack32(WORLD_Y)..pack32(#list)
    for i=1,#list do
        local c = list[i]
        parts[#parts+1] = pack32(c.cx)..pack32(c.cz)
        local bytes = {}
        for j=1,nPer do
            local v = c.data[j] or 0
            if v < 0 then v = 0 elseif v > 255 then v = 255 end
            bytes[#bytes+1] = string.char(v)
        end
        parts[#parts+1] = table.concat(bytes)
    end
    local blob = table.concat(parts)
    gurt.crumbs.set({ name = "world", value = blob })
end

function World:loadFromCrumb()
    local blob = gurt.crumbs.get("world")
    if not blob or #blob < 16 then return false end
    if blob:sub(1,4) ~= "WLD1" then return false end
    local function u32(s, off)
        local b1 = string.byte(s, off) or 0
        local b2 = string.byte(s, off+1) or 0
        local b3 = string.byte(s, off+2) or 0
        local b4 = string.byte(s, off+3) or 0
        local u = b1 + b2*256 + b3*65536 + b4*16777216
        if u >= 2147483648 then u = u - 4294967296 end
        return u, off+4
    end
    local off = 5
    local _, o1 = u32(blob, off); off = o1
    local _, o2 = u32(blob, off); off = o2
    local count; count, off = u32(blob, off)
    local nPer = CHUNK*WORLD_Y*CHUNK
    for i=1,count do
        local cx; cx, off = u32(blob, off)
        local cz; cz, off = u32(blob, off)
        local c = self:getChunk(cx,cz) or self:newChunk(cx,cz)
        for j=1,nPer do
            local pos = off + j - 1
            local val = string.byte(blob, pos) or 0
            c.data[j] = val
        end
        off = off + nPer
    end
    return true
end

local function __hexNibble(n)
    return (n < 10) and (48 + n) or (87 + n)
end

local function __hexByte(n)
    local hi = math.floor(n/16)
    local lo = n % 16
    return string.char(__hexNibble(hi)) .. string.char(__hexNibble(lo))
end

local function __hexToVal(c)
    if c >= 48 and c <= 57 then return c - 48 end
    if c >= 97 and c <= 102 then return c - 87 end
    if c >= 65 and c <= 70 then return c - 55 end
    return 0
end

local function __hexToByte(b1, b2)
    return __hexToVal(b1) * 16 + __hexToVal(b2)
end

local function __hex16(n)
    local hi = math.floor(n / 4096) % 16
    local h2 = math.floor(n / 256) % 16
    local h3 = math.floor(n / 16) % 16
    local lo = n % 16
    return string.char(__hexNibble(hi)) .. string.char(__hexNibble(h2)) .. string.char(__hexNibble(h3)) .. string.char(__hexNibble(lo))
end

local function __hex4ToNum(b1,b2,b3,b4)
    return (__hexToVal(b1) * 4096) + (__hexToVal(b2) * 256) + (__hexToVal(b3) * 16) + __hexToVal(b4)
end

function World:saveChunkToCrumb(cx, cz)
    local c = self:getChunk(cx,cz)
    if not c or not c.dirty then return false end
    local nPer = CHUNK*WORLD_Y*CHUNK
    local t = {"WCR1:"}
    local i = 1
    while i <= nPer do
        local v = c.data[i] or 0
        if v < 0 then v = 0 elseif v > 255 then v = 255 end
        local run = 1
        while i + run <= nPer and c.data[i + run] == v and run < 65535 do
            run = run + 1
        end
        if run >= 3 then
            t[#t+1] = "!" .. __hexByte(v) .. __hex16(run)
            i = i + run
        else
            t[#t+1] = __hexByte(v)
            i = i + 1
        end
    end
    local name = "world:"..cx..","..cz
    gurt.crumbs.set({ name = name, value = table.concat(t) })
    c.dirty = false
    return true
end

function World:loadAllChunksFromCrumbs()
    local all = gurt.crumbs.getAll()
    if not all then return 0 end
    local loaded = 0
    local nPer = CHUNK*WORLD_Y*CHUNK
    for name, crumb in pairs(all) do
        local sx, sz = string.match(name, "^world:(-?%d+),(-?%d+)$")
        if sx and crumb and crumb.value then
            local v = crumb.value
            local cx = tonumber(sx)
            local cz = tonumber(sz)
            local c = self:getChunk(cx,cz) or self:newChunk(cx,cz)
            if #v >= 5 and v:sub(1,5) == "WCR1:" then
                local pos = 6
                local i = 1
                while i <= nPer and pos <= #v do
                    local ch = string.byte(v, pos)
                    if ch == 33 then
                        local b1 = string.byte(v, pos+1)
                        local b2 = string.byte(v, pos+2)
                        local r1 = string.byte(v, pos+3)
                        local r2 = string.byte(v, pos+4)
                        local r3 = string.byte(v, pos+5)
                        local r4 = string.byte(v, pos+6)
                        if not b1 or not b2 or not r1 or not r2 or not r3 or not r4 then break end
                        local val = __hexToByte(b1, b2)
                        local run = __hex4ToNum(r1,r2,r3,r4)
                        for k2=1,run do
                            if i > nPer then break end
                            c.data[i] = val
                            i = i + 1
                        end
                        pos = pos + 7
                    else
                        local b1 = string.byte(v, pos)
                        local b2 = string.byte(v, pos+1)
                        if not b1 or not b2 then break end
                        c.data[i] = __hexToByte(b1, b2)
                        i = i + 1
                        pos = pos + 2
                    end
                end
                c.dirty = false
                loaded = loaded + 1
            end
        end
    end
    return loaded
end

function World:meshChunk(chunkX, chunkZ)
    if not self:getChunk(chunkX - 1, chunkZ) then self:generateChunk(chunkX - 1, chunkZ) end
    if not self:getChunk(chunkX + 1, chunkZ) then self:generateChunk(chunkX + 1, chunkZ) end
    if not self:getChunk(chunkX, chunkZ - 1) then self:generateChunk(chunkX, chunkZ - 1) end
    if not self:getChunk(chunkX, chunkZ + 1) then self:generateChunk(chunkX, chunkZ + 1) end
    local vertices, faces = {}, {}
    local function addVertex(x, y, z)
        vertices[#vertices + 1] = { x * SCALE, y * SCALE, z * SCALE }
        return #vertices
    end
    local height = WORLD_Y
    local cols = CHUNK
    local baseX = chunkX * CHUNK
    local baseZ = chunkZ * CHUNK
    local dirMask, matMask = {}, {}
    local function sample(ix, iy, iz)
        local wx = baseX + (ix - 1)
        local wy = iy - 1
        local wz = baseZ + (iz - 1)
        return self:getBlock(wx, wy, wz)
    end

    for xFace = 0, cols do
        local idx = 1
        local idx = 1
        for y = 1, height do
            for z = 1, cols do
                local aId = sample(xFace, y, z)
                local bId = sample(xFace + 1, y, z)
                local aSolid = aId ~= 0
                local bSolid = bId ~= 0
                if aSolid ~= bSolid then
                    local dir = bSolid and 1 or -1
                    dirMask[idx] = dir
                    matMask[idx] = dir == 1 and bId or aId
                else
                    dirMask[idx] = 0
                    matMask[idx] = 0
                end
                idx = idx + 1
            end
        end
        local rows, widthMax = height, cols
        idx = 1
        for row = 1, rows do
            local col = 1
            while col <= widthMax do
                local dir = dirMask[idx]
                if dir ~= 0 then
                    local mat = matMask[idx]
                    local width = 1
                    while col + width <= widthMax and dirMask[idx + width] == dir and matMask[idx + width] == mat do width = width + 1 end
                    local heightRun = 1
                    while row + heightRun <= rows do
                        local ok = true
                        local rowBase = idx + heightRun * widthMax
                        for c = 0, width - 1 do if dirMask[rowBase + c] ~= dir or matMask[rowBase + c] ~= mat then ok = false break end end
                        if not ok then break end
                        heightRun = heightRun + 1
                    end
                    local x = baseX + xFace
                    local y0 = row - 1
                    local y1 = y0 + heightRun
                    local z0 = baseZ + (col - 1)
                    local z1 = z0 + width
                    local a0 = addVertex(x, y0, z0)
                    local a1 = addVertex(x, y1, z0)
                    local a2 = addVertex(x, y1, z1)
                    local a3 = addVertex(x, y0, z1)
                    local color = MATERIAL_COLOR[mat]
                    if dir == 1 then faces[#faces+1] = {a0,a3,a2,a1,color} else faces[#faces+1] = {a0,a1,a2,a3,color} end
                    for rr = 0, heightRun - 1 do for cc = 0, width - 1 do dirMask[idx + cc + rr * widthMax] = 0; matMask[idx + cc + rr * widthMax] = 0 end end
                    col = col + width
                    idx = idx + width
                else
                    col = col + 1
                    idx = idx + 1
                end
            end
        end
    end

    for yFace = 0, height do
        local idx = 1
        local idx = 1
        for z = 1, cols do
            for x = 1, cols do
                local aId = sample(x, yFace, z)
                local bId = sample(x, yFace + 1, z)
                local aSolid = aId ~= 0
                local bSolid = bId ~= 0
                if aSolid ~= bSolid then
                    local dir = bSolid and 1 or -1
                    dirMask[idx] = dir
                    matMask[idx] = dir == 1 and bId or aId
                else
                    dirMask[idx] = 0
                    matMask[idx] = 0
                end
                idx = idx + 1
            end
        end
        local rows, widthMax = cols, cols
        idx = 1
        for row = 1, rows do
            local col = 1
            while col <= widthMax do
                local dir = dirMask[idx]
                if dir ~= 0 then
                    local mat = matMask[idx]
                    local width = 1
                    while col + width <= widthMax and dirMask[idx + width] == dir and matMask[idx + width] == mat do width = width + 1 end
                    local heightRun = 1
                    while row + heightRun <= rows do
                        local ok = true
                        local rowBase = idx + heightRun * widthMax
                        for c = 0, width - 1 do if dirMask[rowBase + c] ~= dir or matMask[rowBase + c] ~= mat then ok = false break end end
                        if not ok then break end
                        heightRun = heightRun + 1
                    end
                    local y = yFace
                    local x0 = baseX + (col - 1)
                    local x1 = x0 + width
                    local z0 = baseZ + (row - 1)
                    local z1 = z0 + heightRun
                    local a0 = addVertex(x0, y, z0)
                    local a1 = addVertex(x1, y, z0)
                    local a2 = addVertex(x1, y, z1)
                    local a3 = addVertex(x0, y, z1)
                    local color = MATERIAL_COLOR[mat]
                    if dir == 1 then faces[#faces+1] = {a0,a1,a2,a3,color} else faces[#faces+1] = {a0,a3,a2,a1,color} end
                    for rr = 0, heightRun - 1 do for cc = 0, width - 1 do dirMask[idx + cc + rr * widthMax] = 0; matMask[idx + cc + rr * widthMax] = 0 end end
                    col = col + width
                    idx = idx + width
                else
                    col = col + 1
                    idx = idx + 1
                end
            end
        end
    end

    for zFace = 0, cols do
        local idx = 1
        local idx = 1
        for y = 1, height do
            for x = 1, cols do
                local aId = sample(x, y, zFace)
                local bId = sample(x, y, zFace + 1)
                local aSolid = aId ~= 0
                local bSolid = bId ~= 0
                if aSolid ~= bSolid then
                    local dir = bSolid and 1 or -1
                    dirMask[idx] = dir
                    matMask[idx] = dir == 1 and bId or aId
                else
                    dirMask[idx] = 0
                    matMask[idx] = 0
                end
                idx = idx + 1
            end
        end
        local rows, widthMax = height, cols
        idx = 1
        for row = 1, rows do
            local col = 1
            while col <= widthMax do
                local dir = dirMask[idx]
                if dir ~= 0 then
                    local mat = matMask[idx]
                    local width = 1
                    while col + width <= widthMax and dirMask[idx + width] == dir and matMask[idx + width] == mat do width = width + 1 end
                    local heightRun = 1
                    while row + heightRun <= rows do
                        local ok = true
                        local rowBase = idx + heightRun * widthMax
                        for c = 0, width - 1 do if dirMask[rowBase + c] ~= dir or matMask[rowBase + c] ~= mat then ok = false break end end
                        if not ok then break end
                        heightRun = heightRun + 1
                    end
                    local z = baseZ + zFace
                    local x0 = baseX + (col - 1)
                    local x1 = x0 + width
                    local y0 = row - 1
                    local y1 = y0 + heightRun
                    local a0 = addVertex(x0, y0, z)
                    local a1 = addVertex(x1, y0, z)
                    local a2 = addVertex(x1, y1, z)
                    local a3 = addVertex(x0, y1, z)
                    local color = MATERIAL_COLOR[mat]
                    if dir == 1 then faces[#faces+1] = {a0,a3,a2,a1,color} else faces[#faces+1] = {a0,a1,a2,a3,color} end
                    for rr = 0, heightRun - 1 do for cc = 0, width - 1 do dirMask[idx + cc + rr * widthMax] = 0; matMask[idx + cc + rr * widthMax] = 0 end end
                    col = col + width
                    idx = idx + width
                else
                    col = col + 1
                    idx = idx + 1
                end
            end
        end
    end

    local chunk = self:getChunk(chunkX, chunkZ) or self:newChunk(chunkX, chunkZ)
    chunk.verts, chunk.faces = vertices, faces
    return vertices, faces
end

function World:meshChunk2(cx, cz)
    local verts, faces = {}, {}
    local function v(x,y,z)
        verts[#verts+1] = {x*SCALE, y*SCALE, z*SCALE}
        return #verts
    end
    local function quad(x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3)
        local a=v(x0,y0,z0)
        local b=v(x1,y1,z1)
        local c=v(x2,y2,z2)
        local d=v(x3,y3,z3)
        faces[#faces+1]={a,b,c,d}
    end

    for lx=0,CHUNK-1 do
        for ly=0,WORLD_Y-1 do
            for lz=0,CHUNK-1 do
                local wx=cx*CHUNK+lx
                local wy=ly
                local wz=cz*CHUNK+lz
                if self:isSolid(wx,wy,wz) then
                    if not self:isSolid(wx-1,wy,wz) then
                        quad(wx,wy,wz,   wx,wy+1,wz,   wx,wy+1,wz+1,   wx,wy,wz+1)
                    end
                    if not self:isSolid(wx+1,wy,wz) then
                        quad(wx+1,wy,wz,   wx+1,wy,wz+1,   wx+1,wy+1,wz+1,   wx+1,wy+1,wz)
                    end
                    if not self:isSolid(wx,wy-1,wz) then
                        quad(wx,wy,wz,   wx+1,wy,wz,   wx+1,wy,wz+1,   wx,wy,wz+1)
                    end
                    if not self:isSolid(wx,wy+1,wz) then
                        quad(wx,wy+1,wz,   wx,wy+1,wz+1,   wx+1,wy+1,wz+1,   wx+1,wy+1,wz)
                    end
                    if not self:isSolid(wx,wy,wz-1) then
                        quad(wx,wy,wz,   wx,wy+1,wz,   wx+1,wy+1,wz,   wx+1,wy,wz)
                    end
                    if not self:isSolid(wx,wy,wz+1) then
                        quad(wx,wy,wz+1,   wx+1,wy,wz+1,   wx+1,wy+1,wz+1,   wx,wy+1,wz+1)
                    end
                end
            end
        end
    end

    local c = self:getChunk(cx,cz) or self:newChunk(cx,cz)
    c.verts, c.faces = verts, faces
    return verts, faces
end


local RENDER_DISTANCE = 2

function World:updateChunksAroundPlayer(playerX, playerZ)
    local pcx = math.floor(playerX / CHUNK)
    local pcz = math.floor(playerZ / CHUNK)
    local cx0, cz0 = pcx - RENDER_DISTANCE, pcz - RENDER_DISTANCE
    local cx1, cz1 = pcx + RENDER_DISTANCE, pcz + RENDER_DISTANCE

    local wanted = {}
    for cz = cz0, cz1 do
        for cx = cx0, cx1 do
            wanted[k(cx,cz)] = true
        end
    end

    for cz = cz0, cz1 do
        for cx = cx0, cx1 do
            if not self:getChunk(cx,cz) then
                self:generateChunk(cx,cz)
            end
        end
    end

    for cz = cz0, cz1 do
        for cx = cx0, cx1 do
            local verts, faces = self:meshChunk(cx,cz)
            local key = k(cx,cz)
            local mesh = loadedMeshes[key]
            if mesh then
                mesh.vertices = verts
                mesh.faces = faces
            else
                loadedMeshes[key] = Renderer:addMesh(verts, faces, {
                    position = {0,0,0},
                    rotation = {0,0,0},
                    scale = {1,1,1},
                    color = "#00aa00"
                })
            end
        end
    end

    for key, mesh in pairs(loadedMeshes) do
        if not wanted[key] then
            for i, m in ipairs(Renderer.meshes) do
                if m == mesh then
                    table.remove(Renderer.meshes, i)
                    break
                end
            end
            loadedMeshes[key] = nil
        end
    end
end


function World:generateChunk(cx, cz)
    local chunk = self:getChunk(cx, cz)
    if chunk then return end
    
    self:newChunk(cx, cz)
    local bx = cx * CHUNK
    local bz = cz * CHUNK
    
    for lx = 0, CHUNK - 1 do
        for lz = 0, CHUNK - 1 do
            local wx = bx + lx
            local wz = bz + lz
            
            local heightNoise = math.noise(wx * 0.02, wz * 0.02, 0)
            local detailNoise = math.noise(wx * 0.1, wz * 0.1, 0) * 0.3
            local ridgeNoise = math.abs(math.noise(wx * 0.05, wz * 0.05, 100)) * 0.4
            
            local combinedHeight = heightNoise + detailNoise + ridgeNoise
            combinedHeight = (combinedHeight + 1) * 0.5
            
            local h = math.floor(combinedHeight * (WORLD_Y - 1))
            for y = 0, h do
                if y == h then
                    self:setBlock(wx, y, wz, 1)
                else
                    self:setBlock(wx, y, wz, 2)
                end
            end
        end
    end
end


World:loadAllChunksFromCrumbs()
World:updateChunksAroundPlayer(Renderer.camera.pos[1], Renderer.camera.pos[3])
local framequeued = false
local lastPlayerCX, lastPlayerCZ = nil, nil
local function frame()
    if framequeued then framequeued = false onNextFrame(frame) return end
    framequeued = true
    local moved = false
    if keys['W'] then
        local forward = { -math.sin(Renderer.camera.rot[2]), 0, math.cos(Renderer.camera.rot[2]) }
        Renderer.camera.pos[1] = Renderer.camera.pos[1] + forward[1] * 0.2
        Renderer.camera.pos[3] = Renderer.camera.pos[3] + forward[3] * 0.2
        moved = true
    end
    if keys['S'] then
        local forward = { -math.sin(Renderer.camera.rot[2]), 0, math.cos(Renderer.camera.rot[2]) }
        Renderer.camera.pos[1] = Renderer.camera.pos[1] - forward[1] * 0.2
        Renderer.camera.pos[3] = Renderer.camera.pos[3] - forward[3] * 0.2
        moved = true
    end
    if keys['A'] then
        local right = { math.cos(Renderer.camera.rot[2]), 0, math.sin(Renderer.camera.rot[2]) }
        Renderer.camera.pos[1] = Renderer.camera.pos[1] - right[1] * 0.2
        Renderer.camera.pos[3] = Renderer.camera.pos[3] - right[3] * 0.2
        moved = true
    end
    if keys['D'] then
        local right = { math.cos(Renderer.camera.rot[2]), 0, math.sin(Renderer.camera.rot[2]) }
        Renderer.camera.pos[1] = Renderer.camera.pos[1] + right[1] * 0.2
        Renderer.camera.pos[3] = Renderer.camera.pos[3] + right[3] * 0.2
        moved = true
    end
    if keys['Q'] then Renderer.camera.pos[2] = Renderer.camera.pos[2] + 0.2 moved = true end
    if keys['E'] then Renderer.camera.pos[2] = Renderer.camera.pos[2] - 0.2 moved = true end
    if keys['Left'] then Renderer.camera.rot[2] = Renderer.camera.rot[2] + 0.05 moved = true end
    if keys['Right'] then Renderer.camera.rot[2] = Renderer.camera.rot[2] - 0.05 moved = true end
    if keys['Up'] then Renderer.camera.rot[1] = Renderer.camera.rot[1] - 0.05 moved = true end
    if keys['Down'] then Renderer.camera.rot[1] = Renderer.camera.rot[1] + 0.05 moved = true end

    if keys['F'] and Renderer.hitFace then
        local hit = Renderer.hitFace
        World:removeBlock(hit.x, hit.y, hit.z)
        keys['F'] = false
        moved = true
    end
    
    if keys['R'] and Renderer.hitFace then
        local hit = Renderer.hitFace
        local placeX = hit.x + hit.nx
        local placeY = hit.y + hit.ny
        local placeZ = hit.z + hit.nz
        World:placeBlock(placeX, placeY, placeZ, currentBlockId)
        keys['R'] = false
        moved = true
    end
    if keys['1'] then currentBlockId = 1 keys['1'] = false end
    if keys['2'] then currentBlockId = 2 keys['2'] = false end
    
    local currentPlayerCX = math.floor(Renderer.camera.pos[1] / CHUNK)
    local currentPlayerCZ = math.floor(Renderer.camera.pos[3] / CHUNK)
    
    if lastPlayerCX ~= currentPlayerCX or lastPlayerCZ ~= currentPlayerCZ then
        World:updateChunksAroundPlayer(Renderer.camera.pos[1], Renderer.camera.pos[3])
        lastPlayerCX = currentPlayerCX
        lastPlayerCZ = currentPlayerCZ
        moved = true
    end
    Renderer:raycast()
    if moved then
        Renderer:render()
    end
    onNextFrame(frame)
end

gurt.body:on('keydown', function(event) keys[event.key] = true end)
gurt.body:on('keyup', function(event) keys[event.key] = false end)
Time.sleep(2.0)
Renderer:raycast()
Renderer:render()
onNextFrame(frame)

