FROM mambaorg/micromamba:2.4.0 AS base

USER root

ARG ENVIRONMENT
ARG PLUGIN_NAME
ARG QIIME2_USER=qiime2
ARG QIIME2_UID=1000
ARG QIIME2_GID=1000

# Keep the user expected by the micromamba entrypoint, but give it the name and
# numeric IDs used by the development container.
RUN groupmod --new-name "${QIIME2_USER}" --gid "${QIIME2_GID}" "${MAMBA_USER}" \
    && usermod --login "${QIIME2_USER}" \
        --uid "${QIIME2_UID}" \
        --gid "${QIIME2_GID}" \
        --home "/home/${QIIME2_USER}" \
        --move-home "${MAMBA_USER}" \
    && sed -i "s/^${MAMBA_USER}$/${QIIME2_USER}/" /etc/arg_mamba_user \
    && chown -R "${QIIME2_UID}:${QIIME2_GID}" "/home/${QIIME2_USER}"

ENV PLUGIN_NAME=$PLUGIN_NAME \
    ENV_NAME=$PLUGIN_NAME \
    MAMBA_USER=$QIIME2_USER \
    MAMBA_USER_ID=$QIIME2_UID \
    MAMBA_USER_GID=$QIIME2_GID
ENV PATH=/opt/conda/envs/${PLUGIN_NAME}/bin:$PATH \
    LC_ALL=C.UTF-8 LANG=C.UTF-8 \
    MPLBACKEND=agg \
    UNIFRAC_USE_GPU=N \
    HOME=/home/${QIIME2_USER} \
    XDG_CONFIG_HOME=/home/${QIIME2_USER}

WORKDIR /home/${QIIME2_USER}
COPY environment.yml .

RUN apt-get update && apt-get install -y --no-install-recommends wget procps make git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN micromamba create --yes -n ${PLUGIN_NAME} --file environment.yml \
    && micromamba clean --all --yes \
    && chown -R "${QIIME2_UID}:${QIIME2_GID}" /opt/conda \
    && chmod -R a+rwX /opt/conda

USER ${QIIME2_USER}

COPY --chown=${QIIME2_UID}:${QIIME2_GID} . ./plugin
RUN micromamba run -n ${PLUGIN_NAME} pip install --no-cache-dir ./plugin

ENV CONDA_PREFIX=/opt/conda/envs/${PLUGIN_NAME}/
RUN micromamba run -n ${PLUGIN_NAME} qiime dev refresh-cache
RUN echo 'eval "$(micromamba shell hook --shell bash)"' >> $HOME/.bashrc \
    && echo "micromamba activate ${PLUGIN_NAME}" >> $HOME/.bashrc
RUN echo "source tab-qiime" >> $HOME/.bashrc


FROM base AS test

LABEL quay.expires-after=4w

RUN micromamba run -n ${PLUGIN_NAME} pip install --no-cache-dir pytest pytest-cov coverage parameterized pytest-xdist
CMD micromamba run -n ${PLUGIN_NAME} make -f ./plugin/Makefile test-cov

FROM base AS prod

# Important: let any UID modify these directories so that
# `docker run -u UID:GID` works
USER root
RUN rm -rf ./plugin \
    && chmod -R a+rwX /home/${QIIME2_USER}

USER ${QIIME2_USER}
