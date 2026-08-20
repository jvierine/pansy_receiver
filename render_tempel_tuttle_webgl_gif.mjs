#!/usr/bin/env node

import { createServer } from "node:http";
import { mkdtemp, readFile, rm } from "node:fs/promises";
import { spawnSync } from "node:child_process";
import { tmpdir } from "node:os";
import { extname, join, resolve } from "node:path";
import { chromium } from "playwright";

const ROOT = resolve(import.meta.dirname);
const CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome";
const PORT = 8766;
const CAPTURE_FPS = 12;
const FRAME_COUNT = 120;
const MODELS = {
  jupiter: "tempel_tuttle_jupiter_only_webgl.gif",
  "jupiter-pr": "tempel_tuttle_jupiter_pr_drag_webgl.gif",
};
const MIME_TYPES = {
  ".html": "text/html; charset=utf-8",
  ".h5": "application/x-hdf5",
  ".js": "text/javascript; charset=utf-8",
};

function requestedModels() {
  const modelIndex = process.argv.indexOf("--model");
  const model = modelIndex >= 0 ? process.argv[modelIndex + 1] : "both";
  if (model === "both") return Object.keys(MODELS);
  if (!(model in MODELS)) {
    throw new Error(`Unknown model ${model}; use jupiter, jupiter-pr, or both`);
  }
  return [model];
}

function run(command, args) {
  const result = spawnSync(command, args, { stdio: "inherit" });
  if (result.status !== 0) throw new Error(`${command} failed with status ${result.status}`);
}

const server = createServer(async (request, response) => {
  try {
    const pathname = new URL(request.url, `http://127.0.0.1:${PORT}`).pathname;
    const relativePath = pathname === "/" ? "55ptempeltuttle_webgl.html" : pathname.slice(1);
    const path = resolve(ROOT, relativePath);
    if (!path.startsWith(`${ROOT}/`)) throw new Error("Invalid path");
    const content = await readFile(path);
    response.writeHead(200, {
      "Content-Type": MIME_TYPES[extname(path)] || "application/octet-stream",
      "Cache-Control": "no-store",
    });
    response.end(content);
  } catch {
    response.writeHead(404);
    response.end("Not found");
  }
});

await new Promise((resolvePromise) => server.listen(PORT, "127.0.0.1", resolvePromise));
const browser = await chromium.launch({
  headless: true,
  executablePath: CHROME,
  args: ["--use-angle=swiftshader", "--enable-webgl"],
});

try {
  for (const model of requestedModels()) {
    const frameDirectory = await mkdtemp(join(tmpdir(), `tempel-tuttle-${model}-`));
    const page = await browser.newPage({ viewport: { width: 960, height: 960 }, deviceScaleFactor: 1 });
    const url = `http://127.0.0.1:${PORT}/55ptempeltuttle_webgl.html?capture=1&model=${model}`;
    await page.goto(url, { waitUntil: "load" });
    await page.waitForFunction(() => window.orbitViewerReady === true, null, { timeout: 30000 });
    const metadata = await page.evaluate(
      (selectedModel) => window.setOrbitCaptureFrame(selectedModel, 0),
      model,
    );

    for (let outputFrame = 0; outputFrame < FRAME_COUNT; outputFrame += 1) {
      const dataFrame = Math.round(outputFrame * (metadata.frameCount - 1) / (FRAME_COUNT - 1));
      await page.evaluate(
        ([selectedModel, selectedFrame]) => window.setOrbitCaptureFrame(selectedModel, selectedFrame),
        [model, dataFrame],
      );
      await page.screenshot({
        path: join(frameDirectory, `frame_${String(outputFrame).padStart(4, "0")}.png`),
      });
    }
    await page.close();

    const palette = join(frameDirectory, "palette.png");
    const framePattern = join(frameDirectory, "frame_%04d.png");
    const output = resolve(ROOT, MODELS[model]);
    run("ffmpeg", [
      "-y", "-framerate", String(CAPTURE_FPS), "-i", framePattern,
      "-vf", "palettegen=max_colors=192:stats_mode=diff", palette,
    ]);
    run("ffmpeg", [
      "-y", "-framerate", String(CAPTURE_FPS), "-i", framePattern,
      "-i", palette,
      "-lavfi", "paletteuse=dither=sierra2_4a:diff_mode=rectangle",
      "-loop", "0", output,
    ]);
    await rm(frameDirectory, { recursive: true, force: true });
    console.log(`Wrote ${output}`);
  }
} finally {
  await browser.close();
  await new Promise((resolvePromise) => server.close(resolvePromise));
}
