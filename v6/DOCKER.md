# CNVision Docker Guide

This guide documents two valid workflows:

1. **Container-first workflow**:
   `run -> clone -> install -> run -> cp -> commit -> push`
2. **Dockerfile build workflow**:
   `git pull -> docker build -> docker push`

All version references use `v#` as a placeholder.

---

## A) Terminal roles

Use two terminals:

- **Terminal 1 (host terminal)**
  - Docker commands: `run`, `cp`, `ps`, `commit`, `push`, `build`, `rm`
  - Optional local Git commands
- **Terminal 2 (inside container shell)**
  - Linux commands: `apt-get`, `git`, `pip`, `python3`

Important: Docker CLI commands run on host terminal only.

---

## B) Container-first workflow (primary path)

### 1) Clean old containers (Terminal 1)

```bash
docker rm -f cnvision-build
```

Optional full cleanup:

```bash
docker rm -f $(docker ps -aq)
```

### 2) Start base container and enter shell (Terminal 1)

```bash
docker run -it -p 8080:5000 --name cnvision-build python:3.11 /bin/bash
```

### 3) Inside container: install git, clone repo, enter version dir (Terminal 2)

```bash
apt-get update && apt-get install -y git
cd /
git clone https://github.com/mliao2710/CNVision
cd /CNVision/v#
```

### 4) Inside container: install dependencies (Terminal 2)

```bash
pip install flask requests
```

Optional extras:

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

## C) Copy large data file from host

`gene2refseq` may be excluded from GitHub due to size.

Terminal 1:

```bash
docker cp ~/Desktop/cnvision/data/gene2refseq cnvision-build:/CNVision/v#/data/gene2refseq
```

Verify running container:

```bash
docker ps
```

---

## D) Freeze and publish image from running container

### 6) Commit container to a version tag (Terminal 1)

```bash
docker commit <container_id> mliao27/cnvision:v#
```

Example:

```bash
docker commit 256b1c16a204 mliao27/cnvision:v#
```

### 7) Push to Docker Hub (Terminal 1)

```bash
docker push mliao27/cnvision:v#
```

---

## E) Verify pushed image

```bash
docker run --rm -it -p 8080:5000 mliao27/cnvision:v# python3 web_app.py
```

Then open:
- `http://localhost:8080`

---

## F) Dockerfile build workflow (alternative, reproducible)

This workflow builds from the current local checkout.

### What "build from local files" means

`docker build ... .` uses the current directory (`.`) as build context.
The resulting image reflects local files at build time (subject to
`.dockerignore`).

It does **not** auto-update from GitHub. A `git pull` step is still required.

### Steps (Terminal 1)

```bash
cd ~/Desktop/cnvision
git pull origin main
docker build -t mliao27/cnvision:v# .
docker push mliao27/cnvision:v#
```

Run for verification:

```bash
docker run --rm -it -p 8080:5000 -v $(pwd)/data:/app/data mliao27/cnvision:v#
```

---

## G) Updating after small Git changes

Two update patterns are valid.

### Pattern A: update inside running container (matches observed flow)

Terminal 2 (inside `/CNVision/v#`):

```bash
git pull
python3 web_app.py
```

Terminal 1 (host):

```bash
docker cp ~/Desktop/cnvision/data/gene2refseq cnvision-build:/CNVision/v#/data/gene2refseq
docker ps
docker commit <container_id> mliao27/cnvision:v#
docker push mliao27/cnvision:v#
```

### Pattern B: update local repo then rebuild with Dockerfile

Terminal 1:

```bash
cd ~/Desktop/cnvision
git config --global credential.helper osxkeychain
git pull origin main
# apply edits if needed
git add .
git commit -m "v#: <small change summary>"
git push origin main
docker build -t mliao27/cnvision:v# .
docker push mliao27/cnvision:v#
```

GitHub auth note for Pattern B:
- GitHub password auth is not supported for Git operations over HTTPS.
- When prompted for password, paste a GitHub Personal Access Token (PAT).
- Username remains the GitHub username.
- The macOS keychain helper stores the token after first successful push.

Optional `latest` tag update:

```bash
docker tag mliao27/cnvision:v# mliao27/cnvision:latest
docker push mliao27/cnvision:latest
```

---

## H) Common errors and fixes

### `fatal: 'origin' does not appear to be a git repository`
Cause: local repo has no remote configured.

Host terminal fix:

```bash
cd ~/Desktop/cnvision
git remote -v
git remote add origin https://github.com/mliao2710/CNVision.git
git remote -v
git pull origin main
```

### `error: No such remote 'origin'` when running `git remote set-url`
Cause: `origin` does not exist yet.

Use `add` first:

```bash
git remote add origin https://github.com/mliao2710/CNVision.git
```

Then `set-url` is available for future changes:

```bash
git remote set-url origin https://github.com/mliao2710/CNVision.git
```

### `git pull` still fails after remote setup
Check that commands are run in the correct repo folder:

```bash
pwd
git rev-parse --show-toplevel
```

Both should point to the CNVision repo path.

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
Wrong directory. Use:

```bash
cd /CNVision/v#
python3 web_app.py
```

### `docker: command not found` inside container
Docker commands must run on host terminal.

### Host port conflict
Use a different host port:

```bash
docker run --rm -it -p 8081:5000 mliao27/cnvision:v# python3 web_app.py
```

Then open `http://localhost:8081`.
