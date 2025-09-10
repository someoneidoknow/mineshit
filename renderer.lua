--!optimize 2
--!native
local canvas = gurt.select("#poly-test")
local ctx = canvas:withContext("2d")
local Renderer = { 
    meshes = {}, 
    camera = { pos = { 0, -5, -3 }, rot = { 0, 0, 0 }, fov = 200, cx = 400, cy = 300 },
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
local CHUNK = 6
local SCREEN_W, SCREEN_H = 800, 600
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

local function drawTri(p1, p2, p3, color)
    ctx:beginPath()
    ctx:moveTo(p1[1], p1[2])
    ctx:lineTo(p2[1], p2[2])
    ctx:lineTo(p3[1], p3[2])
    ctx:closePath()
    ctx:setFillStyle(color)
    ctx:fill()
    --ctx:setStrokeStyle("#000000")
    --ctx:stroke()
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


local WORLD_Y = 8

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
            local a, b, c = vv[f[1]], vv[f[2]], vv[f[3]]
            local wa, wb, wc = wv[f[1]], wv[f[2]], wv[f[3]]
            local v1x, v1y, v1z = wb[1] - wa[1], wb[2] - wa[2], wb[3] - wa[3]
            local v2x, v2y, v2z = wc[1] - wa[1], wc[2] - wa[2], wc[3] - wa[3]
            local nx, ny, nz = cross(v1x, v1y, v1z, v2x, v2y, v2z)
            nx, ny, nz = normalize(nx, ny, nz)
            local intensity = calculateLighting(nx, ny, nz, self.sun)
            local baseColor = m.color or "#888888"
            local litColor = applyLighting(baseColor, intensity)

            local poly = { {a[1],a[2],a[3]}, {b[1],b[2],b[3]}, {c[1],c[2],c[3]} }
            poly = clipFrustumView(poly, self.camera)
            if #poly >= 3 then
                local sp = {}
                for i=1,#poly do
                    local px, py, pz = projectView(self.camera, poly[i][1], poly[i][2], poly[i][3])
                    sp[i] = { px, py, pz }
                end
                for i=2,#sp-1 do
                    local p1, p2, p3 = sp[1], sp[i], sp[i+1]
                    local key = (p1[3] + p2[3] + p3[3]) * (1/3)
                    all[#all+1] = { p1, p2, p3, litColor, key }
                end
            end
        end
    end
    table.sort(all, function(a,b) return a[5] > b[5] end)
    for i=1,#all do
        local t = all[i]
        drawTri(t[1], t[2], t[3], t[4])
    end
    
    if self.hitFace then
    local h = self.hitFace
    local q
    if h.nx ~= 0 then
        local X = h.nx > 0 and h.x + 1 or h.x
        q = { {X,h.y,h.z}, {X,h.y+1,h.z}, {X,h.y+1,h.z+1}, {X,h.y,h.z+1} }
    elseif h.ny ~= 0 then
        local Y = h.ny > 0 and h.y + 1 or h.y
        q = { {h.x,Y,h.z}, {h.x+1,Y,h.z}, {h.x+1,Y,h.z+1}, {h.x,Y,h.z+1} }
    else
        local Z = h.nz > 0 and h.z + 1 or h.z
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
        for i=2,#sp-1 do
            drawTri(sp[1], sp[i], sp[i+1], "#ffff00")
        end
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
    {1,2,3},{1,3,4},
    {5,6,7},{5,7,8},
    {1,5,6},{1,6,2},
    {2,6,7},{2,7,3},
    {3,7,8},{3,8,4},
    {4,8,5},{4,5,1}
}

local cube = Renderer:addMesh(cubeVerts, cubeFaces, { position = { 0, 2, 0 }, rotation = { 0, 0, 0 }, scale = { 1, 1, 1 }, color = "#ff0000" })

local SCALE = 1

local World = {chunks = {}}
local function k(cx,cz) return cx.."|"..cz end
local function idx(x,y,z) return y*CHUNK*CHUNK + z*CHUNK + x + 1 end

function World:getChunk(cx,cz) return self.chunks[k(cx,cz)] end
function World:newChunk(cx,cz)
    local c = {cx=cx, cz=cz, data={}, verts=nil, faces=nil}
    local n = CHUNK*WORLD_Y*CHUNK
    for i=1,n do c.data[i]=0 end
    self.chunks[k(cx,cz)] = c
    Renderer.worldChunks[k(cx,cz)] = c
    return c
end
function World:placeBlock(wx, wy, wz)
    if wy < 0 or wy >= WORLD_Y then return false end
    local cx = math.floor(wx / CHUNK)
    local cz = math.floor(wz / CHUNK)
    local chunk = self:getChunk(cx, cz)
    if not chunk then return false end
    local lx = wx - cx * CHUNK
    local lz = wz - cz * CHUNK
    chunk.data[idx(lx, wy, lz)] = 1
    self:remeshChunk(cx, cz)
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
    self:remeshChunk(cx, cz)
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

function World:meshChunk(cx,cz)
    local verts, faces = {}, {}
    local function v(x,y,z) verts[#verts+1] = {x*SCALE, y*SCALE, z*SCALE}; return #verts end
    local function add4(a,b,c,d) faces[#faces+1]={a,b,c}; faces[#faces+1]={a,c,d} end
    local Y = WORLD_Y
    local mask = {}
    local function S(ix,iy,iz)
        local wx = cx*CHUNK + (ix-1)
        local wy = iy-1
        local wz = cz*CHUNK + (iz-1)
        return self:isSolid(wx,wy,wz)
    end
    for i=0,CHUNK do
        local n=1
        for j=1,Y do
            for k2=1,CHUNK do
                local a = (i>0) and S(i,j,k2) or false
                local b = (i<CHUNK) and S(i+1,j,k2) or false
                mask[n]=(a~=b) and (b and 1 or -1) or 0
                n=n+1
            end
        end
        local nv,nw=Y,CHUNK
        n=1
        for j=1,nv do
            local k1=1
            while k1<=nw do
                local cval=mask[n]
                if cval~=0 then
                    local w=1
                    while k1+w<=nw and w<4 and mask[n+w]==cval do w=w+1 end
                    local h=1
                    local done=false
                    while j+h<=nv and h<4 do
                        for p=0,w-1 do if mask[n+p+h*nw]~=cval then done=true break end end
                        if done then break end
                        h=h+1
                    end
                    local x = cx*CHUNK + i
                    local y0=(j-1)
                    local y1=y0+h
                    local z0=cz*CHUNK + (k1-1)
                    local z1=z0+w
                    local a0=v(x,y0,z0)
                    local a1=v(x,y1,z0)
                    local a2=v(x,y1,z1)
                    local a3=v(x,y0,z1)
                    if cval==1 then add4(a0,a1,a2,a3) else add4(a0,a3,a2,a1) end
                    for q=0,h-1 do for p=0,w-1 do mask[n+p+q*nw]=0 end end
                    k1=k1+w
                    n=n+w
                else
                    k1=k1+1
                    n=n+1
                end
            end
        end
    end
    for j=0,Y do
        local n=1
        for k2=1,CHUNK do
            for i2=1,CHUNK do
                local a=(j>0) and S(i2,j,k2) or false
                local b=(j<Y) and S(i2,j+1,k2) or false
                mask[n]=(a~=b) and (b and 1 or -1) or 0
                n=n+1
            end
        end
        local nv,nw=CHUNK,CHUNK
        n=1
        for k1=1,nv do
            local i1=1
            while i1<=nw do
                local cval=mask[n]
                if cval~=0 then
                    local w=1
                    while i1+w<=nw and w<4 and mask[n+w]==cval do w=w+1 end
                    local h=1
                    local done=false
                    while k1+h<=nv and h<4 do
                        for p=0,w-1 do if mask[n+p+h*nw]~=cval then done=true break end end
                        if done then break end
                        h=h+1
                    end
                    local y=j
                    local x0=cx*CHUNK + (i1-1)
                    local x1=x0+w
                    local z0=cz*CHUNK + (k1-1)
                    local z1=z0+h
                    local a0=v(x0,y,z0)
                    local a1=v(x1,y,z0)
                    local a2=v(x1,y,z1)
                    local a3=v(x0,y,z1)
                    if cval==1 then add4(a0,a1,a2,a3) else add4(a0,a3,a2,a1) end
                    for q=0,h-1 do for p=0,w-1 do mask[n+p+q*nw]=0 end end
                    i1=i1+w
                    n=n+w
                else
                    i1=i1+1
                    n=n+1
                end
            end
        end
    end
    for k=0,CHUNK do
        local n=1
        for j=1,Y do
            for i2=1,CHUNK do
                local a=(k>0) and S(i2,j,k) or false
                local b=(k<CHUNK) and S(i2,j,k+1) or false
                mask[n]=(a~=b) and (b and 1 or -1) or 0
                n=n+1
            end
        end
        local nv,nw=Y,CHUNK
        n=1
        for j1=1,nv do
            local i1=1
            while i1<=nw do
                local cval=mask[n]
                if cval~=0 then
                    local w=1
                    while i1+w<=nw and w<1 and mask[n+w]==cval do w=w+1 end
                    local h=1
                    local done=false
                    while j1+h<=nv and h<1 do
                        for p=0,w-1 do if mask[n+p+h*nw]~=cval then done=true break end end
                        if done then break end
                        h=h+1
                    end
                    local z=cz*CHUNK + k
                    local x0=cx*CHUNK + (i1-1)
                    local x1=x0+w
                    local y0=(j1-1)
                    local y1=y0+h
                    local a0=v(x0,y0,z)
                    local a1=v(x1,y0,z)
                    local a2=v(x1,y1,z)
                    local a3=v(x0,y1,z)
                    if cval==1 then add4(a0,a1,a2,a3) else add4(a0,a3,a2,a1) end
                    for q=0,h-1 do for p=0,w-1 do mask[n+p+q*nw]=0 end end
                    i1=i1+w
                    n=n+w
                else
                    i1=i1+1
                    n=n+1
                end
            end
        end
    end
    local c = self:getChunk(cx,cz) or self:newChunk(cx,cz)
    c.verts, c.faces = verts, faces
    return verts, faces
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
        faces[#faces+1]={a,b,c}
        faces[#faces+1]={a,c,d}
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


local RENDER_DISTANCE = 1

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
                self:setBlock(wx, y, wz, 1)
            end
        end
    end
end

World:updateChunksAroundPlayer(Renderer.camera.pos[1], Renderer.camera.pos[3])
local framequeued = false
local lastPlayerCX, lastPlayerCZ = nil, nil

local function frame()
    if framequeued then return end
    framequeued = true
    if keys['W'] then
        local forward = { -math.sin(Renderer.camera.rot[2]), 0, math.cos(Renderer.camera.rot[2]) }
        Renderer.camera.pos[1] = Renderer.camera.pos[1] + forward[1] * 0.2
        Renderer.camera.pos[3] = Renderer.camera.pos[3] + forward[3] * 0.2
    end
    if keys['S'] then
        local forward = { -math.sin(Renderer.camera.rot[2]), 0, math.cos(Renderer.camera.rot[2]) }
        Renderer.camera.pos[1] = Renderer.camera.pos[1] - forward[1] * 0.2
        Renderer.camera.pos[3] = Renderer.camera.pos[3] - forward[3] * 0.2
    end
    if keys['A'] then
        local right = { math.cos(Renderer.camera.rot[2]), 0, math.sin(Renderer.camera.rot[2]) }
        Renderer.camera.pos[1] = Renderer.camera.pos[1] - right[1] * 0.2
        Renderer.camera.pos[3] = Renderer.camera.pos[3] - right[3] * 0.2
    end
    if keys['D'] then
        local right = { math.cos(Renderer.camera.rot[2]), 0, math.sin(Renderer.camera.rot[2]) }
        Renderer.camera.pos[1] = Renderer.camera.pos[1] + right[1] * 0.2
        Renderer.camera.pos[3] = Renderer.camera.pos[3] + right[3] * 0.2
    end
    if keys['Q'] then Renderer.camera.pos[2] = Renderer.camera.pos[2] + 0.2 end
    if keys['E'] then Renderer.camera.pos[2] = Renderer.camera.pos[2] - 0.2 end
    if keys['Left'] then Renderer.camera.rot[2] = Renderer.camera.rot[2] + 0.05 end
    if keys['Right'] then Renderer.camera.rot[2] = Renderer.camera.rot[2] - 0.05 end
    if keys['Up'] then Renderer.camera.rot[1] = Renderer.camera.rot[1] - 0.05 end
    if keys['Down'] then Renderer.camera.rot[1] = Renderer.camera.rot[1] + 0.05 end
    
    if keys['F'] and Renderer.hitFace then
        local hit = Renderer.hitFace
        World:removeBlock(hit.x, hit.y, hit.z)
        keys['F'] = false
    end
    
    if keys['R'] and Renderer.hitFace then
        local hit = Renderer.hitFace
        local placeX = hit.x + hit.nx
        local placeY = hit.y + hit.ny
        local placeZ = hit.z + hit.nz
        World:placeBlock(placeX, placeY, placeZ)
        keys['R'] = false
    end
    
    local currentPlayerCX = math.floor(Renderer.camera.pos[1] / CHUNK)
    local currentPlayerCZ = math.floor(Renderer.camera.pos[3] / CHUNK)
    
    if lastPlayerCX ~= currentPlayerCX or lastPlayerCZ ~= currentPlayerCZ then
        World:updateChunksAroundPlayer(Renderer.camera.pos[1], Renderer.camera.pos[3])
        lastPlayerCX = currentPlayerCX
        lastPlayerCZ = currentPlayerCZ
    end
    
    Renderer:raycast()
    Renderer:render()
    framequeued = false
    onNextFrame(frame)
end

gurt.body:on('keydown', function(event) keys[event.key] = true end)
gurt.body:on('keyup', function(event) keys[event.key] = false end)
Time.sleep(2.0)
onNextFrame(frame)


