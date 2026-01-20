# Building and sharing a portable image for cnvision

Quick commands to build an image locally (on your Mac):

Build the image locally (default platform of the Docker engine):

```bash
docker build -t cnvision:latest .
```

# create builder once (if not already present)
docker buildx create --use --name mybuilder
docker buildx build --platform linux/amd64 -t cnvision:latest --load .
```

Share the image by saving it to a tar and sending the file (or via USB/cloud):

```bash
docker save cnvision:latest -o cnvision_latest.tar
docker load -i cnvision_latest.tar
docker run -it --rm -p 5000:5000 cnvision:latest
```

Or push to a registry (Docker Hub example):

```bash
docker tag cnvision:latest mliao27/cnvision:latest
docker push mliao27/cnvision:latest
```

Running with a mounted local data folder (so data isn't baked into the image):

```bash
docker run -it --rm -v $(pwd)/data:/app/data -p 5000:5000 cnvision:latest
```

Notes:
- The `Dockerfile` uses Python 3.11 (your `pyproject.toml` requires >=3.10). Pin a different Python version in the `FROM` line if needed.
- If your project needs system libraries, add them in the `apt-get install` step of the `Dockerfile`.
- If packaging or installing fails on the image, run the build locally and inspect logs, or add necessary build dependencies.
