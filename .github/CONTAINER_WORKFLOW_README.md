# Container Automation Workflow

This document describes the automated workflow for building and publishing Docker containers for mgnify-pipelines-toolkit.

## Overview

The container automation uses a single workflow that:
1. Builds the Python package from source
2. Creates a Docker container with the locally-built package
3. Pushes the container to Quay.io

The container files (`Dockerfile` and `env.yaml`) are stored in the `container/` directory, allowing local builds without GitHub Actions.

## Workflow: Build and Push Container

**Location**: `.github/workflows/build-push-container.yml`

**Triggers**:
- **Automatic**: When a new release is published
- **Manual**: Workflow dispatch for testing

**Actions**:
1. Checks out the repository code
2. Sets up Python 3.11
3. Builds the Python package using `python -m build`
4. Copies the built wheel file to the `container/` directory
5. Extracts version from `pyproject.toml` (or uses manual override)
6. Sets up Docker Buildx
7. Authenticates with Quay.io
8. Builds the Docker image from `container/Dockerfile`
9. Pushes image to Quay.io with tags:
   - `quay.io/microbiome-informatics/mgnify-pipelines-toolkit:{version}`
   - `quay.io/microbiome-informatics/mgnify-pipelines-toolkit:latest`
10. Cleans up wheel files from the container directory

### Required Secrets

This workflow requires the following GitHub secrets:

- **`QUAY_USERNAME`**: Username or robot account for Quay.io authentication
- **`QUAY_PASSWORD`**: Password or robot token for Quay.io with push permissions

### Quay.io Setup

1. Create a robot account in your Quay.io organization (recommended) or use a service account
2. Grant the robot account write permissions to the `microbiome-informatics/mgnify-pipelines-toolkit` repository
3. Generate a token for the robot account
4. Add the credentials as secrets in this repository:
   - Go to Settings → Secrets and variables → Actions → Secrets
   - Add `QUAY_USERNAME` and `QUAY_PASSWORD`

## Container Files

The container is defined by two files in the `container/` directory:

### `container/Dockerfile`

Uses `mambaorg/micromamba:2.5.0` as the base image and:
- Installs conda dependencies from `env.yaml`
- Installs the locally-built wheel file
- Sets up the environment for production use

### `container/env.yaml`

Defines conda dependencies:
- `pyfastx=2.2.0`
- `python=3.11`
- `pip`

The Python package itself is installed from the locally-built wheel, not from PyPI.

## Building Locally

To build the container on your local machine:

```bash
# Build the Python package
python -m build

# Copy the wheel to the container directory
cp dist/*.whl container/

# Build the Docker image
cd container
docker build -t mgnify-pipelines-toolkit:local .

# Clean up
rm *.whl

# Test the container
docker run -it mgnify-pipelines-toolkit:local bash
```

Inside the container, you can verify the installation:

```bash
python --version  # Should be 3.11
pip list | grep mgnify-pipelines-toolkit
get_mpt_version
```

## Manual Workflow Dispatch

The workflow can be manually triggered for testing or creating custom builds:

### Steps:

1. Navigate to Actions → "Build and Push Container to Quay.io"
2. Click "Run workflow"
3. Select the branch to run from
4. Optional inputs:
   - **`tag`**: Override the container tag from `pyproject.toml` (e.g., `dev`, `test`)
     - Leave empty to use the version from `pyproject.toml`
     - Useful for creating test builds with custom tags
   - **`push`**: Whether to push to Quay.io (default: `true`)
     - Set to `false` for dry-run builds (build only, no push)
     - Useful for testing the build process without publishing
5. Click "Run workflow"

### Use Cases:

- **Dry-run testing**: Set `push: false` to test the build without publishing
- **Development builds**: Use custom `tag` like `dev` or `test-feature-x`
- **Pre-release testing**: Verify the build before creating an official release
- **Rebuilding**: Re-build and push a container without creating a new release

## Complete Workflow Flow

### Automated (on Release):

1. **Developer creates and publishes a release** in GitHub
2. **Workflow automatically triggers**
3. **Python package is built** from the released code
4. **Container is built** with the locally-built package
5. **Container is pushed** to Quay.io with version tag and `latest` tag

### Manual (for Testing):

1. **Developer triggers workflow** from Actions tab
2. **Optionally customize** tag and push behavior
3. **Workflow builds** the package and container
4. **Container is pushed** (if `push: true`) or build-only (if `push: false`)

## Troubleshooting

### Python Build Fails
- Check that `pyproject.toml` is valid
- Ensure all source files are committed
- Review the build step logs for missing dependencies

### Container Build Fails
- Verify `container/Dockerfile` and `container/env.yaml` syntax
- Check that the wheel file was copied successfully
- Ensure conda dependencies are available in configured channels
- Test the build locally to reproduce the issue

### Container Push Fails
- Verify `QUAY_USERNAME` and `QUAY_PASSWORD` secrets are correctly set
- Check that the robot account has write permissions to the Quay.io repository
- Ensure the repository `quay.io/microbiome-informatics/mgnify-pipelines-toolkit` exists
- Check Quay.io's status page for outages

### Version Detection Issues
- Verify `pyproject.toml` has the format: `version = "X.Y.Z"`
- Check that the version follows semantic versioning
- Review the workflow logs for parsing errors

## Testing

Before creating a release, test the container build:

1. **Test locally first**:
   ```bash
   python -m build
   cp dist/*.whl container/
   cd container
   docker build -t test-mgnify-pipelines-toolkit .
   docker run -it test-mgnify-pipelines-toolkit bash
   ```

2. **Test with workflow dispatch** (dry-run):
   - Go to Actions → "Build and Push Container to Quay.io"
   - Set `tag` to `test`
   - Set `push` to `false`
   - Verify the build completes successfully

3. **Test with workflow dispatch** (push to Quay.io):
   - Use `tag` like `test` or `dev`
   - Set `push` to `true`
   - Verify the container is pushed and can be pulled

4. **Verify functionality** in the container:
   ```bash
   docker run -it quay.io/microbiome-informatics/mgnify-pipelines-toolkit:test bash
   # Test commands
   python --version
   get_mpt_version
   # Test some tools
   ```

## Maintenance

### Updating Container Dependencies

To update conda dependencies or the base image:

1. Edit `container/env.yaml` or `container/Dockerfile`
2. Test locally using the local build instructions
3. Commit the changes
4. Test with manual workflow dispatch before creating a release

### Updating the Workflow

When modifying `.github/workflows/build-push-container.yml`:

1. Make changes on a feature branch
2. Test with workflow dispatch from your branch
3. Create a PR and verify the changes work as expected
4. Merge to main

## Security Best Practices

1. **Use Robot Accounts**: Create dedicated robot/service accounts for automation rather than personal accounts
2. **Minimal Permissions**: Grant only the minimum required permissions to tokens
3. **Secret Rotation**: Regularly rotate authentication tokens
4. **Scan Containers**: Consider adding container vulnerability scanning to the workflow
5. **Pin Dependencies**: Use specific versions in `env.yaml` for reproducibility
6. **Audit Logging**: Monitor the Actions logs regularly for any unusual activity
