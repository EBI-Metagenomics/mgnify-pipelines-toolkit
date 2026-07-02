FROM mambaorg/micromamba:2.5.0

COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml

RUN micromamba install -y -n base -f /tmp/env.yaml \
    && micromamba install -y -n base conda-forge::procps-ng \
    && micromamba env export --name base --explicit > environment.lock \
    && echo ">> CONDA_LOCK_START" \
    && cat environment.lock \
    && echo "<< CONDA_LOCK_END" \
    && micromamba clean -a -y

USER root

ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
