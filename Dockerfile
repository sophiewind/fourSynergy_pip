FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="4f7b3e1273324aada14e645b2897b60c68a936145abbd38a667a1fa82480e9c4"

# Step 2: Retrieve conda environments
RUN mkdir -p /conda-envs/563e2d751268aa0a2e85531d5194dea1
COPY envs/align.yml /conda-envs/563e2d751268aa0a2e85531d5194dea1/environment.yaml

RUN mkdir -p /conda-envs/947d75ba2fbee9b718729b349c4a952e
COPY envs/basic4cseq.yml /conda-envs/947d75ba2fbee9b718729b349c4a952e/environment.yaml

RUN mkdir -p /conda-envs/dc59948fcdd38ef4ab023c892fa2d4e5
COPY envs/ensemble.yml /conda-envs/dc59948fcdd38ef4ab023c892fa2d4e5/environment.yaml

RUN mkdir -p /conda-envs/9c7dc458616398f78629d1fe9067f072
COPY envs/fastqc.yml /conda-envs/9c7dc458616398f78629d1fe9067f072/environment.yaml

RUN mkdir -p /conda-envs/47defb3093da4c7ae5c34e0c37e36784
COPY envs/foursig.yml /conda-envs/47defb3093da4c7ae5c34e0c37e36784/environment.yaml

RUN mkdir -p /conda-envs/6c31cd779c6a8254fdf704118e4df818
COPY envs/motif.yml /conda-envs/6c31cd779c6a8254fdf704118e4df818/environment.yaml

RUN mkdir -p /conda-envs/45c1abf63747c33b56c5cdd0f7c73c1b
COPY envs/multiqc.yml /conda-envs/45c1abf63747c33b56c5cdd0f7c73c1b/environment.yaml

RUN mkdir -p /conda-envs/54bfa4aaf46b49934e11ad2bd7ea34c0
COPY envs/peakc.yml /conda-envs/54bfa4aaf46b49934e11ad2bd7ea34c0/environment.yaml

RUN mkdir -p /conda-envs/bc08dca7c2e76d55bb9aaea1e9539bec
COPY envs/r3cseq.yml /conda-envs/bc08dca7c2e76d55bb9aaea1e9539bec/environment.yaml

RUN mkdir -p /conda-envs/a7a38bff03d28b4de3680e63be90b2b9
COPY envs/r4cker.yml /conda-envs/a7a38bff03d28b4de3680e63be90b2b9/environment.yaml

# Step 3: Generate conda environments
RUN conda env create --prefix /conda-envs/563e2d751268aa0a2e85531d5194dea1 --file /conda-envs/563e2d751268aa0a2e85531d5194dea1/environment.yaml && \
    conda env create --prefix /conda-envs/947d75ba2fbee9b718729b349c4a952e --file /conda-envs/947d75ba2fbee9b718729b349c4a952e/environment.yaml && \
    conda env create --prefix /conda-envs/dc59948fcdd38ef4ab023c892fa2d4e5 --file /conda-envs/dc59948fcdd38ef4ab023c892fa2d4e5/environment.yaml && \
    conda env create --prefix /conda-envs/9c7dc458616398f78629d1fe9067f072 --file /conda-envs/9c7dc458616398f78629d1fe9067f072/environment.yaml && \
    conda env create --prefix /conda-envs/47defb3093da4c7ae5c34e0c37e36784 --file /conda-envs/47defb3093da4c7ae5c34e0c37e36784/environment.yaml && \
    conda env create --prefix /conda-envs/6c31cd779c6a8254fdf704118e4df818 --file /conda-envs/6c31cd779c6a8254fdf704118e4df818/environment.yaml && \
    conda env create --prefix /conda-envs/45c1abf63747c33b56c5cdd0f7c73c1b --file /conda-envs/45c1abf63747c33b56c5cdd0f7c73c1b/environment.yaml && \
    conda env create --prefix /conda-envs/54bfa4aaf46b49934e11ad2bd7ea34c0 --file /conda-envs/54bfa4aaf46b49934e11ad2bd7ea34c0/environment.yaml && \
    conda env create --prefix /conda-envs/bc08dca7c2e76d55bb9aaea1e9539bec --file /conda-envs/bc08dca7c2e76d55bb9aaea1e9539bec/environment.yaml && \
    conda env create --prefix /conda-envs/a7a38bff03d28b4de3680e63be90b2b9 --file /conda-envs/a7a38bff03d28b4de3680e63be90b2b9/environment.yaml && \
    conda clean --all -y
