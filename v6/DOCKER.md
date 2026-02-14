# CNVision Docker Guide

This guide documents two valid workflows:

1. **Container-first workflow** :
   `run -> clone -> install -> run -> cp -> commit -> push`
2. **Dockerfile build workflow** :
   `git pull -> docker build -> docker push`

---

## A) Terminal roles

Use two terminals:

- **Terminal 1 (host / Mac terminal)**
  - Docker commands: `run`, `cp`, `ps`, `commit`, `push`, `build`
  - Git commands on your local repo
- **Terminal 2 (inside container shell)**
  - Linux commands, `apt-get`, `pip`, `python3 web_app.py`

Important: Docker CLI commands run on host terminal only.

---

## B) Core v6 path (your current primary workflow)

### 1) Clean old containers (Terminal 1)

```bash
docker rm -f cnvision-build
```

If you want to remove all containers:

```bash
docker rm -f $(docker ps -aq)
```

### 2) Start base container and enter shell (Terminal 1)

```bash
docker run -it -p 8080:5000 --name cnvision-build python:3.11 /bin/bash
```

### 3) Inside container: install git + clone repo (Terminal 2)

```bash
apt-get update && apt-get install -y git
cd /
git clone https://github.com/mliao2710/CNVision
cd /CNVision/v6
```

### 4) Inside container: install Python dependencies (Terminal 2)

```bash
pip install flask requests
```

Optional extras if needed:

```bash
pip install biopython pandas numpy
```

### 5) Inside container: run app (Terminal 2)

```bash
python3 web_app.py
```

Open browser:
- `http://localhost:8080`

---

## C) Copy large data file from host (your `gene2refseq` step)

Because `gene2refseq` may be too large for GitHub, copy it from host.

Terminal 1:

```bash
docker cp ~/Desktop/cnvision/data/gene2refseq cnvision-build:/CNVision/v6/data/gene2refseq
```

Verify container is running:

```bash
docker ps
```

---

## D) Freeze and publish image from running container

### 6) Commit container to Docker image tag (Terminal 1)

```bash
docker commit <container_id> mliao27/cnvision:v6
```

Example:

```bash
docker commit 256b1c16a204 mliao27/cnvision:v6
```

### 7) Push image to Docker Hub (Terminal 1)

```bash
docker push mliao27/cnvision:v6
```

---

## E) Verify pushed image

```bash
docker run --rm -it -p 8080:5000 mliao27/cnvision:v6 python3 web_app.py
```

Then open:
- `http://localhost:8080`

---

## F) Dockerfile build workflow (alternative, reproducible)

This workflow builds directly from your **current local repo files**.

### What “build from local files” means

`docker build ... .` uses the current directory (`.`) as the build context.
So the image reflects exactly what is in your local checkout at build time,
subject to `.dockerignore`.

It does **not** auto-pull from GitHub. You must update local repo first.

### Exact steps (Terminal 1)

1. Go to local repo and update it:

```bash
cd ~/Desktop/cnvision
git pull origin main
```

2. Build image from local files:

```bash
docker build -t mliao27/cnvision:v6 .
```

3. Push built image:

```bash
docker push mliao27/cnvision:v6
```

4. Run and verify:

```bash
docker run --rm -it -p 8080:5000 -v $(pwd)/data:/app/data mliao27/cnvision:v6
```

---

## G) How to update image after small Git changes

Use this every time you make code updates.

### Update cycle (Terminal 1)

```bash
cd ~/Desktop/cnvision
git pull origin main
# edit files
git add .
git commit -m "v6: <small change summary>"
git push origin main
```

Rebuild and republish image:

```bash
docker build -t mliao27/cnvision:v6 .
docker push mliao27/cnvision:v6
```

Optional: also update `latest` tag:

```bash
docker tag mliao27/cnvision:v6 mliao27/cnvision:latest
docker push mliao27/cnvision:latest
```

### Versioning tip

For traceability, use unique tags per release:
- `v6.0.0`, `v6.0.1`, `v6.1.0`

Then keep `latest` as optional convenience.

---

## H) Common errors and exact fixes

### `container name already in use`

```bash
docker rm -f cnvision-build
```

### `ModuleNotFoundError: flask`
Inside container:

```bash
pip install flask requests
```

### `python3: can't open file 'web_app.py'`
You are in wrong directory. Use:

```bash
cd /CNVision/v6
python3 web_app.py
```

### `docker: command not found` inside container
Run Docker commands in host terminal, not container terminal.

### Port already in use on host
Change host port mapping:

```bash
docker run --rm -it -p 8081:5000 mliao27/cnvision:v6 python3 web_app.py
```

Then open `http://localhost:8081`.
