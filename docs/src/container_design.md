# Container Design Notes

This page records the decisions and failure modes that matter when building the
eMCP container image. It is intended as a reference for rebuilding the
`Dockerfile` from scratch or for changing the CASA/container stack later.

## Goals

The container should provide a reproducible eMCP runtime with:

- Modular CASA Python packages.
- WSClean.
- AOFlagger and its Python bindings.
- casacore and casacore measures data.
- Headless plotting with CASA `plotms`.
- Compatibility with Docker, Singularity, and Apptainer.

The main complication is that CASA GUI tools are distributed as AppImages inside
Python wheels, while most pipeline runs are headless and often happen inside
Singularity/Apptainer where runtime FUSE is not reliable.

## Base Operating System

The current image uses `ubuntu:24.04`.

Important consequences:

- Ubuntu 24.04 provides Python 3.12 as the system `python3`.
- Ubuntu 24.04 uses time64 package names for Qt 5 runtime libraries, such as
  `libqt5core5t64`, `libqt5gui5t64`, and `libqt5widgets5t64`.
- If the base image changes, re-check all CASA Python compatibility, package
  names, and binary wheel tags.

Do not assume that CASA versions can be mixed freely across Python versions.
The CASA modular packages are released as a coordinated set.

## CASA Version Selection

CASA modular packages must be selected using the compatibility tables in the
CASA documentation:

- CASA 6.7.2: https://casadocs.readthedocs.io/en/v6.7.2/notebooks/introduction.html
- CASA 6.7.0: https://casadocs.readthedocs.io/en/v6.7.0/notebooks/introduction.html

For `ubuntu:24.04`, the natural choice is CASA 6.7.2 because the base image uses
Python 3.12. The documented CASA 6.7.2 package set includes:

```text
casatools==6.7.2.42
casatasks==6.7.2.42
casashell==6.7.2.42
casatestutils==6.7.2.42
casaconfig==1.4.0
casaplotms==2.7.4
casaviewer==2.4.4
casaplotserver==2.0.3
casatablebrowser==0.0.39
casalogger==1.0.23
casafeather==0.0.27
casampi==0.5.9
```

CASA 6.7.0 is documented for Python 3.10. It may be a good choice for an
Ubuntu/Python 3.10 image, but it is not the right target for the current
Ubuntu 24.04/Python 3.12 image.

## Python Dependency Strategy

The package metadata in `pyproject.toml` remains broad enough for non-container
installations, but broad CASA ranges are risky in the container build because
pip may choose a newer, inconsistent, or temporarily broken upstream wheel.

The container therefore installs the documented CASA package set explicitly
before installing eMCP itself. eMCP is then installed with `--no-deps` so pip
does not re-resolve and replace the CASA stack.

This is intentional:

- The container is a controlled runtime.
- The CASA package set is selected from the CASA docs, not from pip's latest
  resolver result.
- eMCP metadata still describes acceptable package ranges for ordinary installs.

When the explicit CASA set changes, also check that `pyproject.toml` does not
contradict it. For example, CASA 6.7.2 uses `casaconfig==1.4.0`, so the project
metadata must allow that version.

## casaviewer Wheel Hash Issue

During one rebuild, pip selected the `casaviewer==2.4.4` Linux wheel tagged
`manylinux_2_34`, and pip rejected it because the downloaded file did not match
the PyPI hash. The same version also provides a `manylinux_2_28` wheel with a
valid hash.

The Dockerfile therefore installs `casaviewer==2.4.4` using the explicit
`manylinux_2_28` wheel URL with its SHA256 hash. This keeps the documented CASA
6.7.2 version while avoiding the bad upstream artifact.

If PyPI is fixed later, this direct URL can be revisited. Until then, leaving
the version unconstrained can make CI fail again.

## Headless Plotting

`plotms(showgui=False)` still needs an X display. Without one, CASA raises:

```text
ERROR: DISPLAY environment variable is not set! Cannot run plotms.
```

The image includes `xvfb`, and eMCP starts `Xvfb` automatically when `DISPLAY`
is unset. This allows plot generation in headless Docker, Singularity, and
Apptainer runs without changing user commands.

Relevant runtime packages include:

```text
xvfb
libqt5core5t64
libqt5gui5t64
libqt5widgets5t64
libxcb-cursor0
libxcb-xinerama0
libxkbcommon-x11-0
```

The exact package names are Ubuntu-version-specific. Re-check them if the base
image changes.

## CASA AppImages and FUSE

The CASA GUI wheels include AppImages:

```text
casaplotms-x86_64.AppImage
casaplotserver-x86_64.AppImage
casaviewer-x86_64.AppImage
```

Running those AppImages inside Singularity/Apptainer can fail with:

```text
dlopen(): error loading libfuse.so.2
AppImages require FUSE to run.
```

Relying on runtime FUSE is fragile, especially on shared systems. The container
build should extract all installed CASA AppImages at build time with
`--appimage-extract`, remove the original AppImage files, and patch installed
Python launchers to call `squashfs-root/AppRun`.

This avoids requiring FUSE at runtime. The build should fail if no CASA AppImage
files are found, because otherwise the image would silently fall back to the
broken runtime path.

## WSClean, AOFlagger, and casacore

The current Dockerfile builds native tools in a separate builder stage:

```text
casacore v3.8.0
WSClean v3.7
AOFlagger v3.5.0
```

The runtime stage copies `/usr/local` and `/usr/share/casacore` from the builder.
This keeps compilers and source trees out of the final image while preserving
the installed native tools and data.

AOFlagger is built with `ENABLE_GUI=OFF`. That is appropriate for the pipeline
runtime, where AOFlagger is used as a command-line/tool dependency rather than
as a GUI application.

If the native tool versions are changed, verify:

- `wsclean --version`
- `aoflagger --version`
- `python3 -c "import aoflagger; print(aoflagger.__file__)"`
- CASA table access to the casacore data path.

## CASA Data

CASA data is not baked into the image. The image writes a CASA site config:

```python
measurespath = "/root/.casa/data"
datapath = ["/root/.casa/data", "/usr/share/casacore/data"]
data_auto_update = False
measures_auto_update = False
nologger = True
nogui = True
```

The README examples bind host CASA data into `/root/.casa/data`.

This avoids large or changing CASA data downloads during image builds, but it
means users must bind a valid CASA data directory when needed.

## Singularity and Apptainer Runtime Considerations

Users may hit local runtime issues unrelated to the image contents:

- If `squashfuse` is missing, Apptainer may extract the SIF to a temporary
  sandbox.
- If the default temporary filesystem is small, extraction can fail with
  `No space left on device`.
- Users can work around that by setting `APPTAINER_TMPDIR` and
  `APPTAINER_CACHEDIR` to a large filesystem.

This is a host configuration issue, not a problem that can be fixed inside the
container image.

## Verification Checklist

Before pushing a Dockerfile change that triggers the container rebuild, check as
much as possible locally:

```bash
git diff --check
python3 -m compileall -q src/eMCP/plots/eMCP_plots.py
```

Check that the CASA package set resolves for Python 3.12 and Linux wheel tags:

```bash
python3 -m pip install --dry-run --ignore-installed --only-binary=:all: \
  --python-version 312 --implementation cp --abi cp312 \
  --platform manylinux_2_28_x86_64 \
  --platform manylinux2014_x86_64 \
  --platform manylinux_2_17_x86_64 \
  protobuf==3.20 \
  casaconfig==1.4.0 \
  casatools==6.7.2.42 \
  casatasks==6.7.2.42 \
  casaplotms==2.7.4 \
  'casaviewer @ https://files.pythonhosted.org/packages/23/00/d997d5cfb0b8458a4119e8e19483b9015c847a3a8145bda306314c102785/casaviewer-2.4.4-py3-none-manylinux_2_28_x86_64.whl#sha256=f1245490638ca053fa8ad74d589aa616725cbe20d4c740342770222f8aac1efd' \
  casashell==6.7.2.42 \
  casaplotserver==2.0.3 \
  casatestutils==6.7.2.42 \
  casatablebrowser==0.0.39 \
  casalogger==1.0.23 \
  casafeather==0.0.27 \
  casampi==0.5.9
```

If Docker is available, build locally before pushing:

```bash
docker build -t emcp-test .
```

Then test the key runtime paths:

```bash
docker run --rm emcp-test python3 -c "import eMCP; import casatasks; import casatools"
docker run --rm emcp-test python3 -c "import aoflagger; print(aoflagger.__file__)"
docker run --rm emcp-test sh -lc 'find /usr/local -name "*AppImage" -print'
docker run --rm emcp-test sh -lc 'find /usr/local -path "*squashfs-root/AppRun" -print'
docker run --rm emcp-test emcp -h
```

The `find /usr/local -name "*AppImage"` command should print nothing after a
successful build, because all CASA AppImages should have been extracted.

For an end-to-end runtime check, run a small pipeline plotting step inside
Docker or Apptainer with no `DISPLAY` set and confirm that:

- eMCP starts `Xvfb`.
- no FUSE/AppImage error appears.
- plot files are written under `weblog/plots/`.

## Known Limits

The current setup is conservative but not fully future-proof:

- It depends on CASA 6.7.2 package availability and the direct `casaviewer`
  wheel URL remaining accessible.
- It assumes Ubuntu 24.04/Python 3.12. A Python 3.10 image should use a different
  CASA package set.
- It avoids runtime FUSE by extracting AppImages, which increases image size
  compared with leaving compressed AppImages in place.
- It does not solve host Apptainer/Singularity problems such as missing
  `squashfuse` or too-small temporary filesystems.

When changing the base OS, Python version, or CASA version, revisit this document
and re-run the full verification checklist.
