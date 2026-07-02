# Container Automation Workflows

This document describes the automated workflows for creating and publishing Docker containers for mgnify-pipelines-toolkit.

## Overview

The automation consists of two workflows in this repository:

1. **Create Container Update PR**: Creates a PR with updated container files when a new release is published
2. **Build and Push Container**: Builds and pushes the Docker image to Quay.io when the PR is merged

## Workflow 1: Create Container Update PR

**Location**: `.github/workflows/create-container-pr.yml`

**Trigger**: When a new release is published

**Actions**:
1. Extracts the version from `pyproject.toml`
2. Creates/updates `container/` directory with:
   - `Dockerfile` - Container definition using micromamba
   - `env.yaml` - Conda environment with mgnify-pipelines-toolkit and dependencies
3. Creates a pull request to this repository with the updated container files

### Required Secrets

No additional secrets required - uses the default `GITHUB_TOKEN` for PR creation.

## Workflow 2: Build and Push Container

**Location**: `.github/workflows/build-push-container.yml`

**Triggers**:
- New release is published
- Manual workflow dispatch

**Actions**:
1. Extracts version from `pyproject.toml`
2. Sets up Docker Buildx
3. Authenticates with Quay.io
4. Builds the Docker image from `container/Dockerfile`
5. Pushes image to Quay.io with tags:
   - `quay.io/microbiome-informatics/mgnify-pipelines-toolkit:{version}`
   - `quay.io/microbiome-informatics/mgnify-pipelines-toolkit:latest`

### Required Secrets

This workflow requires the following GitHub secrets to be configured:

- **`QUAY_USERNAME`**: Username or robot account for Quay.io authentication
- **`QUAY_PASSWORD`**: Password or robot token for Quay.io with push permissions

### Quay.io Setup

1. Create a robot account in your Quay.io organization (recommended) or use a service account
2. Grant the robot account write permissions to the `microbiome-informatics/mgnify-pipelines-toolkit` repository
3. Generate a token for the robot account
4. Add the credentials as secrets in this repository:
   - Go to Settings → Secrets and variables → Actions → Secrets
   - Add `QUAY_USERNAME` and `QUAY_PASSWORD`

## Manual Workflow Dispatch

Both workflows can be manually triggered if needed:

1. Navigate to the Actions tab in the respective repository
2. Select the workflow
3. Click "Run workflow"
4. Fill in any required inputs

## Security Best Practices

1. **Use Robot Accounts**: Create dedicated robot/service accounts for automation rather than personal accounts
2. **Minimal Permissions**: Grant only the minimum required permissions to tokens
3. **Secret Rotation**: Regularly rotate authentication tokens
4. **Environment Protection**: Consider using GitHub Environments with approval requirements for production deployments
5. **Audit Logging**: Monitor the Actions logs regularly for any unusual activity

## Complete Workflow Flow

1. **Developer creates a release** in GitHub
2. **Both workflows trigger simultaneously**:
   - **Workflow 1** creates a PR with updated `container/Dockerfile` and `container/env.yaml`
   - **Workflow 2** builds and pushes the container to Quay.io
3. **Review and merge the PR** (this updates the container files in the repository for future builds)

## Troubleshooting

### PR Creation Fails
- Check GitHub Actions logs for the `create-container-pr` workflow
- Verify `pyproject.toml` contains a valid version string
- Ensure the workflow has permissions to create PRs (default `GITHUB_TOKEN` should work)

### Container Push Fails
- Verify `QUAY_USERNAME` and `QUAY_PASSWORD` secrets are correctly set
- Check that the robot account has write permissions to the Quay.io repository
- Ensure the repository `quay.io/microbiome-informatics/mgnify-pipelines-toolkit` exists
- Review the Dockerfile syntax for errors
- Check Quay.io's status page for outages

### Version Detection Issues
- Verify `pyproject.toml` has the format: `version = "X.Y.Z"`
- Check that the version follows semantic versioning
- Review the workflow logs for parsing errors

## Testing

Before relying on the automation:

1. **Test PR creation**: Use the manual workflow dispatch to test the PR creation workflow
2. **Review generated files**: Manually verify the generated `container/Dockerfile` and `container/env.yaml`
3. **Test locally**: Build the container locally before merging:
   ```bash
   cd container
   docker build -t test-mgnify-pipelines-toolkit .
   docker run -it test-mgnify-pipelines-toolkit bash
   ```
4. **Verify functionality**: Test that the tools work correctly in the container

## Maintenance

When updating the base container or dependencies:

1. Update the Dockerfile/env.yaml template in `.github/workflows/create-container-pr.yml`
2. Test with a pre-release or manually trigger the workflow
3. Update this documentation if the process changes
4. Consider versioning strategy for major base image updates
