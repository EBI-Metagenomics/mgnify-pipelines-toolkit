FROM mambaorg/micromamba:2.5.0

# Copy environment file and local wheel
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
COPY --chown=$MAMBA_USER:$MAMBA_USER *.whl /tmp/

# Install dependencies and local package
RUN micromamba install -y -n base -f /tmp/env.yaml \
    && micromamba install -y -n base conda-forge::procps-ng \
    && micromamba run -n base pip install /tmp/*.whl \
    && micromamba env export --name base --explicit > /environment.lock \
    && echo ">> CONDA_LOCK_START" \
    && cat /environment.lock \
    && echo "<< CONDA_LOCK_END" \
    && micromamba clean -a -y

USER root

ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
