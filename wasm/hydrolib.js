/*!
* Subset of raylib + clay + ui renderer for hydrogen_orbital
* I stole a lot of shit from musl, zozlib and stuff from the interwebs.
 */

function make_environment(env, mathImports = {}) {
    return new Proxy(env, {
        get(target, prop, receiver) {
            if (Object.prototype.hasOwnProperty.call(mathImports, prop)) {
                const v = mathImports[prop];
                return typeof v === "function" ? v.bind(mathImports) : v;
            }
            if (env[prop] !== undefined) {
                return env[prop].bind(env);
            }
            return (...args) => {
                throw new Error(`NOT IMPLEMENTED: ${prop} ${args}`);
            };
        }
    });
}

function hydrolibWasmLibm() {
    const M = Math;
    return {
        cos:    (x) => M.cos(x),
        sin:    (x) => M.sin(x),
        tan:    (x) => M.tan(x),
        cosf:   (x) => M.cos(x),
        sinf:   (x) => M.sin(x),
        powf:   (x, y) => M.pow(x, y),
        fmaxf:  (x, y) => M.max(x, y),
        fminf:  (x, y) => M.min(x, y),
        fmax:   (x, y) => M.max(x, y),
        atan2f: (y, x) => M.atan2(y, x),
        roundf: (x) => M.round(x),
    };
}

function wasmMalloc(exports, bytes) {
    const m = exports.malloc || exports._malloc || exports.wasm_malloc;
    if (typeof m !== "function") return 0;
    return m(bytes) >>> 0;
}

function hydrolibStaticBaseHref() {
    const loc    = typeof document !== "undefined" && document.baseURI ? document.baseURI : window.location.href;
    const url    = new URL(loc);
    let pathname = url.pathname;
    if (!pathname.endsWith("/")) {
        const lastSegment = pathname.substring(pathname.lastIndexOf("/") + 1);
        if (lastSegment.includes("."))
            pathname = pathname.substring(0, pathname.lastIndexOf("/") + 1);
        else
            pathname = pathname + "/";
    }
    url.pathname = pathname || "/";
    return url.href;
}

function setU32(buf, ptr, v) {
    new Uint32Array(buf, ptr, 1)[0] = v >>> 0;
}

function setI32(buf, ptr, v) {
    new Int32Array(buf, ptr, 1)[0] = v | 0;
}

const FONT_OFF_BASE        = 0;
const FONT_OFF_GLYPH_COUNT = 4;
const FONT_OFF_TEX_ID      = 12;
const FONT_OFF_RECS        = 32;
const FONT_OFF_GLYPHS      = 36;
const GLYPH_STRIDE         = 36;
const GLYPH_OFF_ADVANCE    = 12;

const CAMERA_PERSPECTIVE = 0;

const CAMERA_NEAR = 0.05;
const CAMERA_FAR  = 4000.0;

const MESH_OFF_VERTEX_COUNT   = 0;
const MESH_OFF_TRIANGLE_COUNT = 4;
const MESH_OFF_VERTICES       = 8;
const MESH_OFF_COLORS         = 28;
const MESH_OFF_INDICES        = 32;

function hlReadF32(mem, off) {
    return new Float32Array(mem, off, 1)[0];
}

function hlReadI32(mem, off) {
    return new Int32Array(mem, off, 1)[0];
}

function hlReadU32(mem, off) {
    return new Uint32Array(mem, off, 1)[0];
}

function hlMatMultiply(left, right) {
    const result = {};
    result.m0  = left.m0*right.m0 + left.m1*right.m4 + left.m2*right.m8 + left.m3*right.m12;
    result.m1  = left.m0*right.m1 + left.m1*right.m5 + left.m2*right.m9 + left.m3*right.m13;
    result.m2  = left.m0*right.m2 + left.m1*right.m6 + left.m2*right.m10 + left.m3*right.m14;
    result.m3  = left.m0*right.m3 + left.m1*right.m7 + left.m2*right.m11 + left.m3*right.m15;
    result.m4  = left.m4*right.m0 + left.m5*right.m4 + left.m6*right.m8 + left.m7*right.m12;
    result.m5  = left.m4*right.m1 + left.m5*right.m5 + left.m6*right.m9 + left.m7*right.m13;
    result.m6  = left.m4*right.m2 + left.m5*right.m6 + left.m6*right.m10 + left.m7*right.m14;
    result.m7  = left.m4*right.m3 + left.m5*right.m7 + left.m6*right.m11 + left.m7*right.m15;
    result.m8  = left.m8*right.m0 + left.m9*right.m4 + left.m10*right.m8 + left.m11*right.m12;
    result.m9  = left.m8*right.m1 + left.m9*right.m5 + left.m10*right.m9 + left.m11*right.m13;
    result.m10 = left.m8*right.m2 + left.m9*right.m6 + left.m10*right.m10 + left.m11*right.m14;
    result.m11 = left.m8*right.m3 + left.m9*right.m7 + left.m10*right.m11 + left.m11*right.m15;
    result.m12 = left.m12*right.m0 + left.m13*right.m4 + left.m14*right.m8 + left.m15*right.m12;
    result.m13 = left.m12*right.m1 + left.m13*right.m5 + left.m14*right.m9 + left.m15*right.m13;
    result.m14 = left.m12*right.m2 + left.m13*right.m6 + left.m14*right.m10 + left.m15*right.m14;
    result.m15 = left.m12*right.m3 + left.m13*right.m7 + left.m14*right.m11 + left.m15*right.m15;
    return result;
}

function hlMatPerspectiveRad(fovyRad, aspect, nearPlane, farPlane) {
    const result = {};
    const top    = nearPlane*Math.tan(fovyRad*0.5);
    const bottom = -top;
    const right  = top*aspect;
    const left   = -right;
    const rl     = right - left;
    const tb     = top - bottom;
    const fn     = farPlane - nearPlane;
    result.m0  = (nearPlane*2.0)/rl;
    result.m1  = 0;
    result.m2  = 0;
    result.m3  = 0;
    result.m4  = 0;
    result.m5  = (nearPlane*2.0)/tb;
    result.m6  = 0;
    result.m7  = 0;
    result.m8  = (right + left)/rl;
    result.m9  = (top + bottom)/tb;
    result.m10 = -(farPlane + nearPlane)/fn;
    result.m11 = -1;
    result.m12 = 0;
    result.m13 = 0;
    result.m14 = -(farPlane*nearPlane*2.0)/fn;
    result.m15 = 0;
    return result;
}

function hlMatOrtho(left, right, bottom, top, nearPlane, farPlane) {
    const rl = right - left;
    const tb = top - bottom;
    const fn = farPlane - nearPlane;
    const result = {};
    result.m0  = 2.0/rl;
    result.m1  = 0;
    result.m2  = 0;
    result.m3  = 0;
    result.m4  = 0;
    result.m5  = 2.0/tb;
    result.m6  = 0;
    result.m7  = 0;
    result.m8  = 0;
    result.m9  = 0;
    result.m10 = -2.0/fn;
    result.m11 = 0;
    result.m12 = -(left + right)/rl;
    result.m13 = -(top + bottom)/tb;
    result.m14 = -(farPlane + nearPlane)/fn;
    result.m15 = 1;
    return result;
}

function hlMatLookAt(eye, target, up) {
    let vz  = { x: eye.x - target.x, y: eye.y - target.y, z: eye.z - target.z };
    let len = Math.sqrt(vz.x*vz.x + vz.y*vz.y + vz.z*vz.z);
    if (len === 0) len = 1;
    let il = 1.0/len;
    vz.x *= il;
    vz.y *= il;
    vz.z *= il;
    let vx = {
        x: up.y*vz.z - up.z*vz.y,
        y: up.z*vz.x - up.x*vz.z,
        z: up.x*vz.y - up.y*vz.x,
    };
    len = Math.sqrt(vx.x*vx.x + vx.y*vx.y + vx.z*vx.z);
    if (len === 0) len = 1;
    il = 1.0/len;
    vx.x *= il;
    vx.y *= il;
    vx.z *= il;
    const vy = {
        x: vz.y*vx.z - vz.z*vx.y,
        y: vz.z*vx.x - vz.x*vx.z,
        z: vz.x*vx.y - vz.y*vx.x,
    };
    const result = {};
    result.m0  = vx.x;
    result.m1  = vy.x;
    result.m2  = vz.x;
    result.m3  = 0;
    result.m4  = vx.y;
    result.m5  = vy.y;
    result.m6  = vz.y;
    result.m7  = 0;
    result.m8  = vx.z;
    result.m9  = vy.z;
    result.m10 = vz.z;
    result.m11 = 0;
    result.m12 = -(vx.x*eye.x + vx.y*eye.y + vx.z*eye.z);
    result.m13 = -(vy.x*eye.x + vy.y*eye.y + vy.z*eye.z);
    result.m14 = -(vz.x*eye.x + vz.y*eye.y + vz.z*eye.z);
    result.m15 = 1;
    return result;
}

function hlMatToFloat32(mat) {
    return new Float32Array([
        mat.m0, mat.m1, mat.m2, mat.m3,
        mat.m4, mat.m5, mat.m6, mat.m7,
        mat.m8, mat.m9, mat.m10, mat.m11,
        mat.m12, mat.m13, mat.m14, mat.m15,
    ]);
}

function hlCompileGlProgram(gl, vsSrc, fsSrc) {
    function compile(type, src) {
        const sh = gl.createShader(type);
        gl.shaderSource(sh, src);
        gl.compileShader(sh);
        if (!gl.getShaderParameter(sh, gl.COMPILE_STATUS)) {
            throw new Error(gl.getShaderInfoLog(sh) || "shader compile fail");
        }
        return sh;
    }
    const vs = compile(gl.VERTEX_SHADER, vsSrc);
    const fs = compile(gl.FRAGMENT_SHADER, fsSrc);
    const prog = gl.createProgram();
    gl.attachShader(prog, vs);
    gl.attachShader(prog, fs);
    gl.linkProgram(prog);
    if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) {
        throw new Error(gl.getProgramInfoLog(prog) || "program link fail");
    }
    gl.deleteShader(vs);
    gl.deleteShader(fs);
    return prog;
}

let iota = 0;
const LOG_ALL     = iota++;
const LOG_TRACE   = iota++;
const LOG_DEBUG   = iota++;
const LOG_INFO    = iota++;
const LOG_WARNING = iota++;
const LOG_ERROR   = iota++;
const LOG_FATAL   = iota++;
const LOG_NONE    = iota++;

class HydrolibJs {
    #FONT_SCALE_MAGIC = 0.65;

    #reset() {
        this.previous                   = undefined;
        this.exports                    = undefined;
        this.canvas                     = undefined;
        this.canvasGl                   = undefined;
        this.ctx                        = undefined;
        this.gl                         = undefined;
        this.glProgTri                  = undefined;
        this.glProgPt                   = undefined;
        this.glLocMvpTri                = undefined;
        this.glLocMvpPt                 = undefined;
        this.glLocPtSize                = undefined;
        this.glVAO                      = undefined;
        this.glVBO                      = undefined;
        this._glReady                   = false;
        this._mvpArr                    = undefined;
        this._batchTri                  = undefined;
        this._batchLine                 = undefined;
        this._batchPt                   = undefined;
        this._inMode3D                  = false;
        this._disposeListeners          = undefined;
        this.dt                         = undefined;
        this.targetFPS                  = 60;
        this.entryFunction              = undefined;
        this.prevPressedKeyState        = new Set();
        this.currentPressedKeyState     = new Set();
        this.currentMouseWheelMoveState = 0;
        this.currentMouseWheelMoveV     = [0, 0];
        this.currentMousePosition       = {x: 0, y: 0};
        this.mouseButtonsDown           = new Set();
        this._prevMouseCanvas           = null;
        this._mouseDelta                = [0, 0];
        this.images                     = [];
        this.quit                       = false;
        this.traceLogLevel              = LOG_FATAL;
        this._fontFaceSeq               = 0;
        this._fontFamilyByPtr           = new Map();
        this._scissorStack              = [];
        this._fontBytesCache           = undefined;
    }

    constructor() {
        this.#reset();
    }

    async #preloadFontBytes() {
        this._fontBytesCache = new Map();
        const canonical = "resources/fonts/Iosevka-Regular.ttc";
        const base      = hydrolibStaticBaseHref();
        const relTries = [canonical, "fonts/Iosevka-Regular.ttc"];
        for (const rel of relTries) {
            try {
                const r = await fetch(new URL(rel, base).href);
                if (!r.ok) continue;
                const buf = new Uint8Array(await r.arrayBuffer());
                if (!buf.byteLength) continue;
                this._fontBytesCache.set(canonical, buf);
                return;
            } catch (_) {}
        }
    }

    stop() {
        this.quit = true;
    }

    startExports({ exports, canvasId }) {
        if (this.exports !== undefined) {
            console.error("The game is already running. Please stop() it first.");
            return;
        }

        const canvas = document.getElementById(canvasId);
        this.canvas = canvas;
        this.ctx = canvas.getContext("2d", { alpha: true, desynchronized: false });
        if (this.ctx === null) {
            throw new Error("Could not create 2d canvas context");
        }

        this.canvasGl = document.getElementById("canvas-gl") || undefined;

        this.exports = exports;

        const dispose = [];
        const listen = (target, type, handler, opts) => {
            target.addEventListener(type, handler, opts);
            dispose.push(() => target.removeEventListener(type, handler, opts));
        };

        const keyDown = (e) => {
            this.currentPressedKeyState.add(glfwKeyMapping[e.code]);
        };
        const keyUp = (e) => {
            this.currentPressedKeyState.delete(glfwKeyMapping[e.code]);
        };
        const wheelMove = (e) => {
            this.currentMouseWheelMoveState = Math.sign(-e.deltaY);
            this.currentMouseWheelMoveV[0]  = Math.sign(e.deltaX)  || 0;
            this.currentMouseWheelMoveV[1]  = Math.sign(-e.deltaY) || 0;
        };
        const mouseMove = (e) => {
            this.currentMousePosition = { x: e.clientX, y: e.clientY };
        };
        const mouseDown = (e) => {
            this.mouseButtonsDown.add(e.button);
        };
        const mouseUp = (e) => {
            this.mouseButtonsDown.delete(e.button);
        };

        listen(window, "keydown", keyDown);
        listen(window, "keyup", keyUp);
        listen(window, "wheel", wheelMove, { passive: true });
        listen(window, "mousemove", mouseMove);
        listen(window, "mousedown", mouseDown);
        listen(window, "mouseup", mouseUp);

        this._disposeListeners = dispose;

        this.exports.main();
        const next = (timestamp) => {
            if (this.quit) {
                this.ctx.clearRect(0, 0, this.ctx.canvas.width, this.ctx.canvas.height);
                if (this.gl) {
                    const gl = this.gl;
                    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
                    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
                }
                for (const off of dispose) {
                    off();
                }
                this.#reset();
                return;
            }
            this.dt = (timestamp - this.previous)/1000.0;
            this.previous = timestamp;
            const rootCanvas = this.canvasGl || this.canvas;
            const rect = rootCanvas.getBoundingClientRect();
            const cx = this.currentMousePosition.x - rect.left;
            const cy = this.currentMousePosition.y - rect.top;
            if (this._prevMouseCanvas) {
                this._mouseDelta[0] = cx - this._prevMouseCanvas.x;
                this._mouseDelta[1] = cy - this._prevMouseCanvas.y;
            } else {
                this._mouseDelta[0] = this._mouseDelta[1] = 0;
            }
            this._prevMouseCanvas = { x: cx, y: cy };
            this.entryFunction();
            this.currentMouseWheelMoveState = 0;
            this.currentMouseWheelMoveV[0] = 0;
            this.currentMouseWheelMoveV[1] = 0;
            window.requestAnimationFrame(next);
        };
        window.requestAnimationFrame((timestamp) => {
            this.previous = timestamp;
            window.requestAnimationFrame(next);
        });
    }

    async start({ wasmPath, canvasId }) {
        await this.#preloadFontBytes();
        const wasmHref = new URL(wasmPath || "app.wasm", hydrolibStaticBaseHref()).href;
        const wasm = await WebAssembly.instantiateStreaming(fetch(wasmHref), {
            env: make_environment(this, hydrolibWasmLibm()),
        });
        this.startExports({
            exports: wasm.instance.exports,
            canvasId,
        });
    }

    InitWindow(width, height, title_ptr) {
        this.ctx.canvas.width = width;
        this.ctx.canvas.height = height;
        if (this.canvasGl) {
            this.canvasGl.width = width;
            this.canvasGl.height = height;
        }
        const buffer = this.exports.memory.buffer;
        document.title = cstr_by_ptr(buffer, title_ptr);
    }

    WindowShouldClose() {
        return false;
    }

    SetConfigFlags(_flags) {
        // Do nothing and win
    }

    CloseWindow() {
        this.quit = true;
    }

    SetTraceLogLevel(level) {
        this.traceLogLevel = level;
    }

    SetTargetFPS(fps) {
        console.log(`The game wants to run at ${fps} FPS, but in Web we gonna just ignore it.`);
        this.targetFPS = fps;
    }

    GetScreenWidth() {
        return this.ctx.canvas.width;
    }

    GetScreenHeight() {
        return this.ctx.canvas.height;
    }

    GetFrameTime() {
        // TODO: This is a stopgap solution to prevent sudden jumps in dt when the user switches to a differen tab.
        // We need a proper handling of Target FPS here.
        return Math.min(this.dt, 1.0/this.targetFPS);
    }

    BeginDrawing() {}

    EndDrawing() {
        this.prevPressedKeyState.clear();
        this.prevPressedKeyState = new Set(this.currentPressedKeyState);
    }

    DrawCircleV(center_ptr, radius, color_ptr) {
        const buffer = this.exports.memory.buffer;
        const [x, y] = new Float32Array(buffer, center_ptr, 2);
        const color  = hydrolibRgbaStyleFromPointer(buffer, color_ptr);
        this.ctx.beginPath();
        this.ctx.arc(x, y, radius, 0, 2*Math.PI, false);
        this.ctx.fillStyle = color;
        this.ctx.fill();
    }

    ClearBackground(color_ptr) {
        const buf = this.exports.memory.buffer;
        const rgba = new Uint8Array(buf, color_ptr, 4);
        const rf = rgba[0] / 255;
        const gf = rgba[1] / 255;
        const bf = rgba[2] / 255;

        if (this.canvasGl) {
            try {
                this.#ensureWebgl();
                const gl = this.gl;
                gl.bindFramebuffer(gl.FRAMEBUFFER, null);
                gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
                gl.clearColor(rf, gf, bf, 1.0);
                gl.clearDepth(1.0);
                gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
            } catch (_) {}
            this.ctx.clearRect(0, 0, this.ctx.canvas.width, this.ctx.canvas.height);
            return;
        }

        this.ctx.fillStyle = getColorFromMemory(buf, color_ptr);
        this.ctx.fillRect(0, 0, this.ctx.canvas.width, this.ctx.canvas.height);
    }

    DrawText(text_ptr, posX, posY, fontSize, color_ptr) {
        const buffer = this.exports.memory.buffer;
        const text   = cstr_by_ptr(buffer, text_ptr);
        const color  = getColorFromMemory(buffer, color_ptr);
        fontSize *= this.#FONT_SCALE_MAGIC;
        this.ctx.fillStyle = color;
        this.ctx.font = `${fontSize}px monospace`;

        const lines = text.split('\n');
        for (var i = 0; i < lines.length; i++) {
            this.ctx.fillText(lines[i], posX, posY + fontSize + (i * fontSize));
        }
    }

    DrawRectangle(posX, posY, width, height, color_ptr) {
        const buffer = this.exports.memory.buffer;
        const color  = getColorFromMemory(buffer, color_ptr);
        this.ctx.fillStyle = color;
        this.ctx.fillRect(posX, posY, width, height);
    }

    DrawRectangleV(position_ptr, size_ptr, color_ptr) {
        const buffer   = this.exports.memory.buffer;
        const color    = getColorFromMemory(buffer, color_ptr);
        const position = new Float32Array(buffer, position_ptr, 2);
        const size     = new Float32Array(buffer, size_ptr, 2);
        this.ctx.fillStyle = color;
        this.ctx.fillRect(position[0], position[1], size[0], size[1]);
    }

    IsKeyPressed(key) {
        return !this.prevPressedKeyState.has(key) && this.currentPressedKeyState.has(key);
    }
    IsKeyDown(key) {
        return this.currentPressedKeyState.has(key);
    }
    GetMouseWheelMove() {
      return this.currentMouseWheelMoveState;
    }

    GetMouseWheelMoveV(out_ptr) {
        const b = this.exports.memory.buffer;
        new Float32Array(b, out_ptr, 2).set(this.currentMouseWheelMoveV);
    }

    GetMouseDelta(out_ptr) {
        const b = this.exports.memory.buffer;
        new Float32Array(b, out_ptr, 2).set(this._mouseDelta);
    }

    IsMouseButtonDown(button) {
        return this.mouseButtonsDown.has(button) ? 1 : 0;
    }
    IsGestureDetected() {
        return false;
    }

    TextFormat(... args){
        // TODO: Implement printf style formatting for TextFormat
        return args[0];
    }

    TraceLog(logLevel, text_ptr, ... args) {
        if (logLevel < this.traceLogLevel) return;
        const buffer = this.exports.memory.buffer;
        const text   = cstr_by_ptr(buffer, text_ptr);
        switch(logLevel) {
        case LOG_ALL:     console.log(`ALL: ${text} ${args}`);     break;
        case LOG_TRACE:   console.log(`TRACE: ${text} ${args}`);   break;
        case LOG_DEBUG:   console.log(`DEBUG: ${text} ${args}`);   break;
        case LOG_INFO:    console.log(`INFO: ${text} ${args}`);    break;
        case LOG_WARNING: console.log(`WARNING: ${text} ${args}`); break;
        case LOG_ERROR:   console.log(`ERROR: ${text} ${args}`);   break;
        case LOG_FATAL:   throw new Error(`FATAL: ${text}`);
        case LOG_NONE:    console.log(`NONE: ${text} ${args}`);    break;
        }
    }

    GetMousePosition(result_ptr) {
        const c      = this.ctx.canvas;
        const bcrect = c.getBoundingClientRect();
        const rw     = Math.max(1, bcrect.width);
        const rh     = Math.max(1, bcrect.height);
        const x      = (this.currentMousePosition.x - bcrect.left)*(c.width/rw);
        const y      = (this.currentMousePosition.y - bcrect.top)*(c.height/rh);
        const buffer = this.exports.memory.buffer;
        new Float32Array(buffer, result_ptr, 2).set([x, y]);
    }

    CheckCollisionPointRec(point_ptr, rec_ptr) {
        const buffer           = this.exports.memory.buffer;
        const [x, y]           = new Float32Array(buffer, point_ptr, 2);
        const [rx, ry, rw, rh] = new Float32Array(buffer, rec_ptr, 4);
        return ((x >= rx) && x <= (rx + rw) && (y >= ry) && y <= (ry + rh));
    }

    Fade(result_ptr, color_ptr, alpha) {
        const buffer       = this.exports.memory.buffer;
        const [r, g, b, _] = new Uint8Array(buffer, color_ptr, 4);
        const newA         = Math.max(0, Math.min(255, 255.0*alpha));
        new Uint8Array(buffer, result_ptr, 4).set([r, g, b, newA]);
    }

    DrawRectangleRec(rec_ptr, color_ptr) {
        const buffer       = this.exports.memory.buffer;
        const [x, y, w, h] = new Float32Array(buffer, rec_ptr, 4);
        const color        = getColorFromMemory(buffer, color_ptr);
        this.ctx.fillStyle = color;
        this.ctx.fillRect(x, y, w, h);
    }

    DrawRectangleLinesEx(rec_ptr, lineThick, color_ptr) {
        const buffer       = this.exports.memory.buffer;
        const [x, y, w, h] = new Float32Array(buffer, rec_ptr, 4);
        const color        = getColorFromMemory(buffer, color_ptr);
        this.ctx.strokeStyle = color;
        this.ctx.lineWidth = lineThick;
        this.ctx.strokeRect(x + lineThick/2, y + lineThick/2, w - lineThick, h - lineThick);
    }

    MeasureText(text_ptr, fontSize) {
        const buffer = this.exports.memory.buffer;
        const text   = cstr_by_ptr(buffer, text_ptr);
        fontSize *= this.#FONT_SCALE_MAGIC;
        this.ctx.font = `${fontSize}px monospace`;
        return this.ctx.measureText(text).width;
    }

    TextSubtext(text_ptr, position, length) {
        const buffer  = this.exports.memory.buffer;
        const text    = cstr_by_ptr(buffer, text_ptr);
        const subtext = text.substring(position, length);

        var bytes = new Uint8Array(buffer, 0, subtext.length+1);
        for(var i = 0; i < subtext.length; i++) {
            bytes[i] = subtext.charCodeAt(i);
        }
        bytes[subtext.length] = 0;

        return bytes;
    }

    LoadTexture(result_ptr, filename_ptr) {
        const buffer   = this.exports.memory.buffer;
        const filename = cstr_by_ptr(buffer, filename_ptr);

        var result = new Uint32Array(buffer, result_ptr, 5)
        var img    = new Image();
        img.src = filename;
        this.images.push(img);

        result[0] = this.images.indexOf(img);
        // TODO: get the true width and height of the image
        result[1] = 256; // width
        result[2] = 256; // height
        result[3] = 1; // mipmaps
        result[4] = 7; // format PIXELFORMAT_UNCOMPRESSED_R8G8B8A8

        return result;
    }

    DrawTexture(texture_ptr, posX, posY, color_ptr) {
        const buffer                               = this.exports.memory.buffer;
        const [id, width, height, mipmaps, format] = new Uint32Array(buffer, texture_ptr, 5);
        // // TODO: implement tinting for DrawTexture
        // const tint = getColorFromMemory(buffer, color_ptr);

        this.ctx.drawImage(this.images[id], posX, posY);
    }

    DrawCircle(centerX, centerY, radius, color_ptr) {
        const buffer = this.exports.memory.buffer;
        const color  = getColorFromMemory(buffer, color_ptr);
        this.ctx.beginPath();
        this.ctx.arc(centerX, centerY, radius, 0, 2*Math.PI, false);
        this.ctx.fillStyle = color;
        this.ctx.fill();
    }

    GetFontDefault(result_ptr) {
        this._hydrolibPrepareFontStruct(result_ptr >>> 0, 48, "monospace");
    }

    LoadFontEx(result_ptr, fileName_ptr, fontSize, _codepoints_ptr, _codepointCount) {
        const mem = this.exports.memory.buffer;
        const rel = cstr_by_ptr(mem, fileName_ptr);
        let css   = "monospace";
        const bytes = this._fontBytesCache && this._fontBytesCache.get(rel);
        if (bytes && bytes.byteLength > 0) {
            try {
                const blob    = new Blob([bytes]);
                const blobUrl = URL.createObjectURL(blob);
                css = `hydrolib_font_${++this._fontFaceSeq}`;
                const ff = new FontFace(css, `url(${blobUrl})`);
                document.fonts.add(ff);
                ff.load().catch(() => {});
            } catch (e) {
                console.warn("HYDROLIB->LoadFontEx:", rel, e);
            }
        } else {
            console.warn("HYDROLIB->LoadFontEx: missing preload for", rel, "(serve site from repo root so ../resources works, or adjust preload paths)");
        }
        this._hydrolibPrepareFontStruct(result_ptr >>> 0, fontSize | 0, css);
    }

    _hydrolibPrepareFontStruct(fontPtr, baseSize, cssFamily) {
        // bro, I'm about to kill myself :)
        const exp         = this.exports;
        const glyphsBytes = 95 * GLYPH_STRIDE;
        const recsBytes   = 95 * 16;
        const gp          = wasmMalloc(exp, glyphsBytes);
        const rp          = wasmMalloc(exp, recsBytes);
        const mem = exp.memory.buffer;
        if (!gp || !rp) {
            console.warn("HYDROLIB: wasm malloc not linked. Clay text metrics need libc/allocator in wasm.");
            return;
        }
        const adv = Math.round(baseSize * 0.52);
        const gw = Math.max(1, adv);
        for (let i = 0; i < 95; i++) {
            const g = gp + i*GLYPH_STRIDE;
            setI32(mem, g, 32 + i);
            setI32(mem, g + 4, 0);
            setI32(mem, g + 8, 0);
            setI32(mem, g + GLYPH_OFF_ADVANCE, adv);
            for (let k = 16; k < GLYPH_STRIDE; k += 4) setI32(mem, g + k, 0);
            const rec = rp + i * 16;
            new Float32Array(mem, rec, 4).set([i*gw*0, 0, gw, baseSize]);
        }
        setI32(mem, fontPtr + FONT_OFF_BASE, baseSize);
        setI32(mem, fontPtr + FONT_OFF_GLYPH_COUNT, 95);
        setI32(mem, fontPtr + 8, 2);
        setU32(mem, fontPtr + FONT_OFF_TEX_ID, 0);
        setI32(mem, fontPtr + 16, 1);
        setI32(mem, fontPtr + 20, 1);
        setI32(mem, fontPtr + 24, 1);
        setI32(mem, fontPtr + 28, 7);
        setU32(mem, fontPtr + FONT_OFF_RECS, rp >>> 0);
        setU32(mem, fontPtr + FONT_OFF_GLYPHS, gp >>> 0);
        this._fontFamilyByPtr.set(fontPtr >>> 0, cssFamily);
    }

    _hydrolibFontCss(fontPtr) {
        return this._fontFamilyByPtr.get(fontPtr >>> 0) || "monospace";
    }

    GenTextureMipmaps() {}

    SetTextureFilter() {}

    MeasureTextEx(result_ptr, font_ptr, text_ptr, fontSize, spacing) {
        const buffer = this.exports.memory.buffer;
        const text   = cstr_by_ptr(buffer, text_ptr);
        const out    = new Float32Array(buffer, result_ptr, 2);
        const family = this._hydrolibFontCss(font_ptr);
        const fs     = fontSize * this.#FONT_SCALE_MAGIC;
        this.ctx.font = `${fs}px ${family}`;
        const w = this.ctx.measureText(text).width + Math.max(0, (text.length - 1)*spacing);
        out[0] = w;
        out[1] = fontSize;
    }

    DrawTextEx(font_ptr, text_ptr, position_ptr, fontSize, spacing, tint_ptr) {
        const buffer       = this.exports.memory.buffer;
        const text         = cstr_by_ptr(buffer, text_ptr);
        const [posX, posY] = new Float32Array(buffer, position_ptr, 2);
        const tint         = getColorFromMemory(buffer, tint_ptr);
        const family       = this._hydrolibFontCss(font_ptr);
        const fs           = fontSize*this.#FONT_SCALE_MAGIC;
        this.ctx.fillStyle = tint;
        this.ctx.font      = `${fs}px ${family}`;
        let x              = posX;
        for (let i = 0; i < text.length; i++) {
            const character = text[i];
            if (character === "\n") continue;
            this.ctx.fillText(character, x, posY + fs);
            x += this.ctx.measureText(character).width + spacing;
        }
    }

    DrawTexturePro(texture_ptr, source_rec_ptr, dest_rec_ptr, origin_ptr, rotation, tint_ptr) {
        const b     = this.exports.memory.buffer;
        const texId = new Uint32Array(b, texture_ptr, 1)[0];
        const img   = this.images[texId];
        if (!img || !img.complete || !img.naturalWidth) return;
        const [sx, sy, sw, sh] = new Float32Array(b, source_rec_ptr, 4);
        const [dx, dy, dw, dh] = new Float32Array(b, dest_rec_ptr, 4);
        const [ox, oy]         = new Float32Array(b, origin_ptr, 2);
        this.ctx.save();
        const tint = getColorFromMemory(b, tint_ptr);
        this.ctx.fillStyle = tint;
        this.ctx.translate(dx + ox, dy + oy);
        if (rotation) this.ctx.rotate(rotation*Math.PI/180);
        this.ctx.translate(-ox, -oy);
        this.ctx.drawImage(img, sx, sy, sw, sh, 0, 0, dw, dh);
        this.ctx.restore();
    }

    BeginScissorMode(x, y, w, h) {
        this.ctx.save();
        const iw = w | 0;
        const ih = h | 0;
        if (iw > 0 && ih > 0) {
            this.ctx.beginPath();
            this.ctx.rect(x, y, iw, ih);
            this.ctx.clip();
        }
    }

    EndScissorMode() {
        this.ctx.restore();
    }

    DrawRectangleRounded(rec_ptr, roundness, _segments, color_ptr) {
        const b              = this.exports.memory.buffer;
        const [x, y, rw, rh] = new Float32Array(b, rec_ptr, 4);
        const col            = getColorFromMemory(b, color_ptr);
        const rr             = Number.isFinite(roundness) ? roundness : 0;
        const r              = Math.max(0, Math.min(rr, rw/2, rh/2));
        this.ctx.save();
        this.ctx.fillStyle = col;
        let drew = false;
        if (typeof this.ctx.roundRect === "function" && r > 0 && rw > 0 && rh > 0) {
            try {
                this.ctx.beginPath();
                this.ctx.roundRect(x, y, rw, rh, r);
                this.ctx.fill();
                drew = true;
            } catch (_) {}
        }
        if (!drew) this.ctx.fillRect(x, y, rw, rh);
        this.ctx.restore();
    }

    DrawRing(center_ptr, innerR, outerR, startDeg, endDeg, _segments, color_ptr) {
        const b        = this.exports.memory.buffer;
        const [cx, cy] = new Float32Array(b, center_ptr, 2);
        const col      = getColorFromMemory(b, color_ptr);
        const mid      = (innerR + outerR)/2;
        const w        = Math.max(1, outerR - innerR);
        this.ctx.save();
        this.ctx.strokeStyle = col;
        this.ctx.lineWidth = w;
        this.ctx.beginPath();
        this.ctx.arc(cx, cy, mid, startDeg*Math.PI/180, endDeg*Math.PI/180);
        this.ctx.stroke();
        this.ctx.restore();
    }

    #ensureWebgl() {
        if (this._glReady || !this.canvasGl) return;
        const gl = this.canvasGl.getContext("webgl2", {
            alpha: false,
            antialias: false,
            depth: true,
            stencil: false,
            preserveDrawingBuffer: false,
        });
        if (!gl) {
            throw new Error("WebGL2 unavailable");
        }
        this.gl = gl;
        const vsTri = `#version 300 es
uniform mat4 uMVP;
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec4 aCol;
out vec4 vCol;
void main() {
    vCol = aCol;
    gl_Position = uMVP * vec4(aPos, 1.0);
}`;
        const fs = `#version 300 es
precision highp float;
in vec4 vCol;
out vec4 fragColor;
void main() { fragColor = vCol; }`;
        const vsPt = `#version 300 es
uniform mat4 uMVP;
uniform float uPointSize;
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec4 aCol;
out vec4 vCol;
void main() {
    vCol = aCol;
    gl_PointSize = uPointSize;
    gl_Position = uMVP * vec4(aPos, 1.0);
}`;
        this.glProgTri = hlCompileGlProgram(gl, vsTri, fs);
        this.glProgPt = hlCompileGlProgram(gl, vsPt, fs);
        this.glLocMvpTri = gl.getUniformLocation(this.glProgTri, "uMVP");
        this.glLocMvpPt = gl.getUniformLocation(this.glProgPt, "uMVP");
        this.glLocPtSize = gl.getUniformLocation(this.glProgPt, "uPointSize");
        this.glVAO = gl.createVertexArray();
        gl.bindVertexArray(this.glVAO);
        this.glVBO = gl.createBuffer();
        gl.bindBuffer(gl.ARRAY_BUFFER, this.glVBO);
        gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 28, 0);
        gl.vertexAttribPointer(1, 4, gl.FLOAT, false, 28, 12);
        gl.enableVertexAttribArray(0);
        gl.enableVertexAttribArray(1);
        gl.bindVertexArray(null);
        gl.enable(gl.DEPTH_TEST);
        gl.depthFunc(gl.LEQUAL);
        this._glReady = true;
    }

    #computeMvp(camera_ptr) {
        const mem = this.exports.memory.buffer;
        const px = hlReadF32(mem, camera_ptr);
        const py = hlReadF32(mem, camera_ptr + 4);
        const pz = hlReadF32(mem, camera_ptr + 8);
        const tx = hlReadF32(mem, camera_ptr + 12);
        const ty = hlReadF32(mem, camera_ptr + 16);
        const tz = hlReadF32(mem, camera_ptr + 20);
        const ux = hlReadF32(mem, camera_ptr + 24);
        const uy = hlReadF32(mem, camera_ptr + 28);
        const uz = hlReadF32(mem, camera_ptr + 32);
        const fovyDeg = hlReadF32(mem, camera_ptr + 36);
        const projection = hlReadI32(mem, camera_ptr + 40);
        const w = this.canvasGl ? this.canvasGl.width : this.ctx.canvas.width;
        const h = this.canvasGl ? this.canvasGl.height : this.ctx.canvas.height;
        const aspect = w / Math.max(h, 1);
        const view = hlMatLookAt({ x: px, y: py, z: pz }, { x: tx, y: ty, z: tz }, { x: ux, y: uy, z: uz });
        let proj;
        if (projection === CAMERA_PERSPECTIVE) {
            proj = hlMatPerspectiveRad(fovyDeg * Math.PI / 180.0, aspect, CAMERA_NEAR, CAMERA_FAR);
        } else {
            const top = fovyDeg / 2.0;
            const right = top * aspect;
            proj = hlMatOrtho(-right, right, -top, top, CAMERA_NEAR, CAMERA_FAR);
        }
        const mvp = hlMatMultiply(view, proj);
        if (!this._mvpArr) this._mvpArr = new Float32Array(16);
        this._mvpArr.set(hlMatToFloat32(mvp));
    }

    #flushGlBatches() {
        const gl = this.gl;
        if (!gl || !this._mvpArr) return;
        gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
        gl.bindVertexArray(this.glVAO);
        gl.bindBuffer(gl.ARRAY_BUFFER, this.glVBO);
        gl.useProgram(this.glProgTri);
        gl.uniformMatrix4fv(this.glLocMvpTri, false, this._mvpArr);
        if (this._batchTri.length > 0) {
            const n = this._batchTri.length / 7;
            gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(this._batchTri), gl.STREAM_DRAW);
            gl.drawArrays(gl.TRIANGLES, 0, n);
        }
        if (this._batchLine.length > 0) {
            const n = this._batchLine.length / 7;
            gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(this._batchLine), gl.STREAM_DRAW);
            gl.drawArrays(gl.LINES, 0, n);
        }
        if (this._batchPt.length > 0) {
            gl.useProgram(this.glProgPt);
            gl.uniformMatrix4fv(this.glLocMvpPt, false, this._mvpArr);
            gl.uniform1f(this.glLocPtSize, 4.0);
            const n = this._batchPt.length / 7;
            gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(this._batchPt), gl.STREAM_DRAW);
            gl.drawArrays(gl.POINTS, 0, n);
        }
        gl.bindVertexArray(null);
    }

    BeginMode3D(camera_ptr) {
        if (!this.canvasGl) return;
        try {
            this.#ensureWebgl();
        } catch (_) {
            return;
        }
        this._batchTri = [];
        this._batchLine = [];
        this._batchPt = [];
        this.#computeMvp(camera_ptr);
        this._inMode3D = true;
    }

    EndMode3D() {
        if (!this.gl || !this._inMode3D) {
            this._inMode3D = false;
            return;
        }
        this._inMode3D = false;
        this.#flushGlBatches();
    }

    DrawModel(_model_ptr, _posX, _posY, _posZ, _scale, _tint_ptr) {}

    DrawLine3D(start_ptr, end_ptr, color_ptr) {
        const b = this.exports.memory.buffer;
        const s = new Float32Array(b, start_ptr, 3);
        const e = new Float32Array(b, end_ptr, 3);
        const rgba = new Uint8Array(b, color_ptr, 4);
        const cr = rgba[0] / 255;
        const cg = rgba[1] / 255;
        const cb = rgba[2] / 255;
        const ca = rgba[3] / 255;
        if (this.gl && this._inMode3D) {
            const a = this._batchLine;
            a.push(s[0], s[1], s[2], cr, cg, cb, ca, e[0], e[1], e[2], cr, cg, cb, ca);
            return;
        }
        const c = getColorFromMemory(b, color_ptr);
        this.ctx.save();
        this.ctx.strokeStyle = c;
        this.ctx.lineWidth = 1;
        this.ctx.beginPath();
        this.ctx.moveTo(s[0], s[1]);
        this.ctx.lineTo(e[0], e[1]);
        this.ctx.stroke();
        this.ctx.restore();
    }

    DrawPoint3D(pos_ptr, color_ptr) {
        const b = this.exports.memory.buffer;
        const p = new Float32Array(b, pos_ptr, 3);
        const rgba = new Uint8Array(b, color_ptr, 4);
        const cr = rgba[0] / 255;
        const cg = rgba[1] / 255;
        const cb = rgba[2] / 255;
        const ca = rgba[3] / 255;
        if (this.gl && this._inMode3D) {
            this._batchPt.push(p[0], p[1], p[2], cr, cg, cb, ca);
            return;
        }
        const c = getColorFromMemory(b, color_ptr);
        this.ctx.fillStyle = c;
        this.ctx.fillRect(p[0], p[1], 2, 2);
    }

    DrawMeshIndexedWeb(mesh_ptr, camera_ptr) {
        if (!this.gl || !this._inMode3D || !mesh_ptr) return;
        const mem = this.exports.memory.buffer;
        const triCount = hlReadI32(mem, mesh_ptr + MESH_OFF_TRIANGLE_COUNT);
        const vc = hlReadI32(mem, mesh_ptr + MESH_OFF_VERTEX_COUNT);
        const vertsPtr = hlReadU32(mem, mesh_ptr + MESH_OFF_VERTICES);
        const colsPtr = hlReadU32(mem, mesh_ptr + MESH_OFF_COLORS);
        const idxPtr = hlReadU32(mem, mesh_ptr + MESH_OFF_INDICES);
        if (triCount <= 0 || vc <= 0 || !vertsPtr || !colsPtr || !idxPtr) return;

        const camx = hlReadF32(mem, camera_ptr);
        const camy = hlReadF32(mem, camera_ptr + 4);
        const camz = hlReadF32(mem, camera_ptr + 8);

        const verts = new Float32Array(mem, vertsPtr, vc * 3);
        const colors = new Uint8Array(mem, colsPtr, vc * 4);
        const indices = new Uint16Array(mem, idxPtr, triCount * 3);
        const tri = this._batchTri;

        for (let i = 0; i < triCount; i++) {
            const i0 = indices[i * 3];
            const i1 = indices[i * 3 + 1];
            const i2 = indices[i * 3 + 2];
            if (i0 >= vc || i1 >= vc || i2 >= vc) continue;

            const v0x = verts[i0 * 3], v0y = verts[i0 * 3 + 1], v0z = verts[i0 * 3 + 2];
            const v1x = verts[i1 * 3], v1y = verts[i1 * 3 + 1], v1z = verts[i1 * 3 + 2];
            const v2x = verts[i2 * 3], v2y = verts[i2 * 3 + 1], v2z = verts[i2 * 3 + 2];

            const e1x = v1x - v0x, e1y = v1y - v0y, e1z = v1z - v0z;
            const e2x = v2x - v0x, e2y = v2y - v0y, e2z = v2z - v0z;
            let nx = e2y * e1z - e2z * e1y;
            let ny = e2z * e1x - e2x * e1z;
            let nz = e2x * e1y - e2y * e1x;
            let nlen = Math.sqrt(nx * nx + ny * ny + nz * nz);
            if (nlen === 0) continue;
            nx /= nlen;
            ny /= nlen;
            nz /= nlen;

            const cx = (v0x + v1x + v2x) / 3;
            const cy = (v0y + v1y + v2y) / 3;
            const cz = (v0z + v1z + v2z) / 3;
            let fx = camx - cx, fy = camy - cy, fz = camz - cz;
            let flen = Math.sqrt(fx * fx + fy * fy + fz * fz);
            if (flen === 0) continue;
            fx /= flen;
            fy /= flen;
            fz /= flen;

            if (nx * fx + ny * fy + nz * fz <= 0) continue;

            const pushV = (vi) => {
                const j = vi * 3;
                const k = vi * 4;
                tri.push(
                    verts[j], verts[j + 1], verts[j + 2],
                    colors[k] / 255, colors[k + 1] / 255, colors[k + 2] / 255, colors[k + 3] / 255,
                );
            };
            pushV(i0);
            pushV(i2);
            pushV(i1);
        }
    }

    LoadMaterialDefault(result_ptr) {
        const b = this.exports.memory.buffer;
        const n = 512;
        new Uint8Array(b, result_ptr, n).fill(0);
    }

    UploadMesh(_mesh_ptr, _dynamic) {}

    UnloadMaterial(_material_ptr) {}

    UnloadMesh(_mesh_ptr) {}

    GetRandomValue(min, max) {
        return min + Math.floor(Math.random()*(max - min + 1));
    }

    IsGamepadAvailable(_gamepad) {
        return false;
    }

    GetGamepadAxisMovement(_gamepad, _axis) {
        return 0;
    }

    ColorFromHSV(result_ptr, hue, saturation, value) {
        const buffer = this.exports.memory.buffer;
        const result = new Uint8Array(buffer, result_ptr, 4);

        let k = (5.0 + hue/60.0)%6;
        let t = 4.0 - k;
        k = (t < k) ? t : k;
        k = (k < 1) ? k : 1;
        k = (k > 0) ? k : 0;
        result[0] = Math.floor((value - value*saturation*k)*255.0);

        k = (3.0 + hue/60.0)%6;
        t = 4.0 - k;
        k = (t < k) ? t : k;
        k = (k < 1) ? k : 1;
        k = (k > 0) ? k : 0;
        result[1] = Math.floor((value - value*saturation*k)*255.0);

        k = (1.0 + hue/60.0)%6;
        t = 4.0 - k;
        k = (t < k) ? t : k;
        k = (k < 1) ? k : 1;
        k = (k > 0) ? k : 0;
        result[2] = Math.floor((value - value*saturation*k)*255.0);

        result[3] = 255;
    }

    hydrolib_js_set_entry(entry) {
        this.entryFunction = this.exports.__indirect_function_table.get(entry);
    }
}

globalThis.HydrolibJs = HydrolibJs;

const glfwKeyMapping = {
    "Space":          32,
    "Quote":          39,
    "Comma":          44,
    "Minus":          45,
    "Period":         46,
    "Slash":          47,
    "Digit0":         48,
    "Digit1":         49,
    "Digit2":         50,
    "Digit3":         51,
    "Digit4":         52,
    "Digit5":         53,
    "Digit6":         54,
    "Digit7":         55,
    "Digit8":         56,
    "Digit9":         57,
    "Semicolon":      59,
    "Equal":          61,
    "KeyA":           65,
    "KeyB":           66,
    "KeyC":           67,
    "KeyD":           68,
    "KeyE":           69,
    "KeyF":           70,
    "KeyG":           71,
    "KeyH":           72,
    "KeyI":           73,
    "KeyJ":           74,
    "KeyK":           75,
    "KeyL":           76,
    "KeyM":           77,
    "KeyN":           78,
    "KeyO":           79,
    "KeyP":           80,
    "KeyQ":           81,
    "KeyR":           82,
    "KeyS":           83,
    "KeyT":           84,
    "KeyU":           85,
    "KeyV":           86,
    "KeyW":           87,
    "KeyX":           88,
    "KeyY":           89,
    "KeyZ":           90,
    "BracketLeft":    91,
    "Backslash":      92,
    "BracketRight":   93,
    "Backquote":      96,
    //  GLFW_KEY_WORLD_1   161
    //  GLFW_KEY_WORLD_2   162
    "Escape":         256,
    "Enter":          257,
    "Tab":            258,
    "Backspace":      259,
    "Insert":         260,
    "Delete":         261,
    "ArrowRight":     262,
    "ArrowLeft":      263,
    "ArrowDown":      264,
    "ArrowUp":        265,
    "PageUp":         266,
    "PageDown":       267,
    "Home":           268,
    "End":            269,
    "CapsLock":       280,
    "ScrollLock":     281,
    "NumLock":        282,
    "PrintScreen":    283,
    "Pause":          284,
    "F1":             290,
    "F2":             291,
    "F3":             292,
    "F4":             293,
    "F5":             294,
    "F6":             295,
    "F7":             296,
    "F8":             297,
    "F9":             298,
    "F10":            299,
    "F11":            300,
    "F12":            301,
    "F13":            302,
    "F14":            303,
    "F15":            304,
    "F16":            305,
    "F17":            306,
    "F18":            307,
    "F19":            308,
    "F20":            309,
    "F21":            310,
    "F22":            311,
    "F23":            312,
    "F24":            313,
    "F25":            314,
    "NumPad0":        320,
    "NumPad1":        321,
    "NumPad2":        322,
    "NumPad3":        323,
    "NumPad4":        324,
    "NumPad5":        325,
    "NumPad6":        326,
    "NumPad7":        327,
    "NumPad8":        328,
    "NumPad9":        329,
    "NumpadDecimal":  330,
    "NumpadDivide":   331,
    "NumpadMultiply": 332,
    "NumpadSubtract": 333,
    "NumpadAdd":      334,
    "NumpadEnter":    335,
    "NumpadEqual":    336,
    "ShiftLeft":      340,
    "ControlLeft" :   341,
    "AltLeft":        342,
    "MetaLeft":       343,
    "ShiftRight":     344,
    "ControlRight":   345,
    "AltRight":       346,
    "MetaRight":      347,
    "ContextMenu":    348,
    //  GLFW_KEY_LAST   GLFW_KEY_MENU
}

function cstrlen(mem, ptr) {
    let len = 0;
    while (mem[ptr] != 0) {
        len++;
        ptr++;
    }
    return len;
}

function cstr_by_ptr(mem_buffer, ptr) {
    const mem   = new Uint8Array(mem_buffer);
    const len   = cstrlen(mem, ptr);
    const bytes = new Uint8Array(mem_buffer, ptr, len);
    return new TextDecoder().decode(bytes);
}

function hydrolibRgbaStyleFromPointer(buffer, color_ptr) {
    const [r, g, b, a] = new Uint8Array(buffer, color_ptr, 4);
    return `rgba(${r},${g},${b},${a/255})`;
}

function getColorFromMemory(buffer, color_ptr) {
    return hydrolibRgbaStyleFromPointer(buffer, color_ptr);
}
