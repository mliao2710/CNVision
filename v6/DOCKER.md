# CNVision Docker + Two-Terminal Runbook

This runbook matches your actual workflow and common errors seen in real usage.
Use it exactly as written when building/running/testing/pushing.

## Terminal roles (important)

Use two terminals:

- **Terminal A (Mac host)**
  - `git clone/pull/add/commit/push`
  - `docker build/run/exec/tag/push`
- **Terminal B (inside container shell)**
  - optional debugging and running `python3 web_app.py` in repo path

Rule: run Docker commands only in Terminal A (host), not inside container.

## 0) Quick sanity checks (Terminal A)

```bash
pwd
ls
```

You should be in your repo root (example: `/Users/milesliao/Desktop/cnvision`).

## 1) Git workflow (Terminal A)

```bash
git pull origin main
git status
# edit files
git add .
git commit -m "v6: <message>"
git push origin main
```

## 2) Build and run your app image (Terminal A)

Preferred: run the image you build locally.

```bash
docker build -t cnvision:latest .
docker run --rm -it -p 8080:5000 -v $(pwd)/data:/app/data cnvision:latest
```

Open:
- `http://localhost:8080`

## 3) If you want a shell inside your app container (Terminal A -> Terminal B)

Start container in background:

```bash
docker run -d --name cnvision-run -p 8080:5000 -v $(pwd)/data:/app/data cnvision:latest
```

Enter shell:

```bash
docker exec -it cnvision-run /bin/bash
```

Inside container (Terminal B):

```bash
pwd
ls
python3 web_app.py
```

## 4) If you use raw python base image (debug style you used)

If you start from `python:3.11`, you must install deps and clone manually.
This is for debugging only, not normal release.

```bash
docker run -it -p 8080:5000 --name cnvision-build python:3.11 /bin/bash
```

Inside container:

```bash
apt-get update && apt-get install -y git
cd /
git clone https://github.com/mliao2710/CNVision
cd /CNVision/v5   # adjust if using v6 folder structure
pip install flask requests
python3 web_app.py
```

## 5) Docker Hub push (Terminal A)

```bash
docker tag cnvision:latest mliao27/cnvision:v6
docker push mliao27/cnvision:v6
```

Optional latest tag:

```bash
docker tag cnvision:latest mliao27/cnvision:latest
docker push mliao27/cnvision:latest
```

## 6) Common errors you hit and exact fixes

### Error: `/bin/ba: no such file or directory`
Cause: typo.

Use:

```bash
/bin/bash
```

### Error: image tag not found (`mliao27/cnvision:v5: not found`)
Cause: tag does not exist remotely yet.

Fix:
- build locally and tag/push first, or
- run an existing tag.

### Error: container name already in use (`/cnvision-build` conflict)
Fix:

```bash
docker rm -f cnvision-build
```

Then rerun.

### Error: `Address already in use` (port 5000 busy on Mac)
Fix option A (recommended): map to a different host port:

```bash
docker run --rm -it -p 8080:5000 cnvision:latest
```

Fix option B: stop process using 5000.

### Error: `docker: command not found` inside container
Cause: running Docker command inside container shell.

Fix:
- run Docker commands in Terminal A (host), not Terminal B.

### Error: `python3: can't open file '//web_app.py'`
Cause: wrong working directory in container.

Fix:

```bash
cd /CNVision/v5
python3 web_app.py
```

## 7) Stop/cleanup commands (Terminal A)

Stop one container:

```bash
docker stop cnvision-run
```

Remove one container:

```bash
docker rm -f cnvision-run
```

Remove many old containers:

```bash
docker container prune -f
```

List containers/images:

```bash
docker ps -a
docker images | grep cnvision
```

## 8) Gold path (copy/paste, minimal friction)

Terminal A:

```bash
cd ~/Desktop/cnvision
git pull origin main
docker build -t cnvision:latest .
docker run --rm -it -p 8080:5000 -v $(pwd)/data:/app/data cnvision:latest
```

Then open `http://localhost:8080`.

If release push:

```bash
git add .
git commit -m "v6: final"
git push origin main
docker tag cnvision:latest mliao27/cnvision:v6
docker push mliao27/cnvision:v6
```
