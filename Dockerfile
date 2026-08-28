FROM continuumio/miniconda3:latest AS base

USER root
ARG ENVIRONMENT
ARG PLUGIN_NAME
ARG QIIME2_USER=qiime2
ARG QIIME2_UID=1000
ARG QIIME2_GID=1000
RUN groupadd --gid "${QIIME2_GID}" "${QIIME2_USER}" \
    && useradd --uid "${QIIME2_UID}" \
        --gid "${QIIME2_GID}" \
        --create-home \
        --shell /bin/bash \
        "${QIIME2_USER}"
ENV PLUGIN_NAME=$PLUGIN_NAME \
    ENV_NAME=$PLUGIN_NAME
ENV PATH=/opt/conda/envs/${PLUGIN_NAME}/bin:$PATH \
    LC_ALL=C.UTF-8 LANG=C.UTF-8 \
    MPLBACKEND=agg \
    UNIFRAC_USE_GPU=N \
    HOME=/home/${QIIME2_USER} \
    XDG_CONFIG_HOME=/home/${QIIME2_USER}

WORKDIR /home/${QIIME2_USER}
COPY environment.yml .

RUN apt-get update && apt-get install -y --no-install-recommends wget procps make \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# The environment is created and assigned to the runtime user in one layer, so
# changing its ownership does not duplicate it in the final image.
RUN conda update -qy conda \
    && conda install -c conda-forge -qy mamba=2.4.0 \
    && mamba env create -n ${PLUGIN_NAME} --file environment.yml \
    && mamba clean --all --yes \
    && chown -R "${QIIME2_UID}:${QIIME2_GID}" "/opt/conda/envs/${PLUGIN_NAME}" \
    && chmod -R a+rwX "/opt/conda/envs/${PLUGIN_NAME}" \

USER ${QIIME2_USER}

COPY --chown=${QIIME2_UID}:${QIIME2_GID} . ./plugin

RUN mamba run -n ${PLUGIN_NAME} pip install --no-cache-dir ./plugin

ENV CONDA_PREFIX=/opt/conda/envs/${PLUGIN_NAME}/
RUN mamba run -n ${PLUGIN_NAME} qiime dev refresh-cache
RUN echo "source activate ${PLUGIN_NAME}" >> $HOME/.bashrc \
    && echo "source tab-qiime" >> $HOME/.bashrc
FROM base AS test

LABEL quay.expires-after=4w
RUN mamba run -n ${PLUGIN_NAME} pip install --no-cache-dir pytest pytest-cov coverage parameterized pytest-xdist
CMD ["/bin/bash", "-c", "exec mamba run -n \"${PLUGIN_NAME}\" make -f ./plugin/Makefile test-cov"]

FROM base AS prod
# Important: let any UID modify this directory so that
# `docker run -u UID:GID` works.
USER root
RUN rm -rf ./plugin \
    && chmod -R a+rwX /home/${QIIME2_USER} \

USER ${QIIME2_USER}
