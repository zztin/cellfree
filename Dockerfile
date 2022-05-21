# syntax=docker/dockerfile:1
FROM python:3.9.13-bullseye AS development_build
LABEL maintainer="litingchen16@gmail.com"

ARG CELLFREE_ENV
ENV CELLFREE_ENV=${CELLFREE_ENV} \
  PYTHONFAULTHANDLER=1 \
  PYTHONUNBUFFERED=1 \
  PYTHONHASHSEED=random \
  PIP_NO_CACHE_DIR=off \
  PIP_DISABLE_PIP_VERSION_CHECK=on \
  PIP_DEFAULT_TIMEOUT=100 \
  POETRY_VERSION=1.1.3
# System development tool
RUN pip install "poetry==$POETRY_VERSION"
# Copy only requirements to cache them in docker layer
WORKDIR /cellfree-workspace
COPY poetry.lock pyproject.toml /cellfree-workspace/
# Project initialization:
RUN poetry config virtualenvs.create false
RUN poetry install --no-root --no-dev --no-interaction --no-ansi
# Creating folders, and files for a project:
COPY ./cellfree-0.0.2-py3-none-any.whl /cellfree-workspace/
RUN pip install /cellfree-workspace/cellfree-0.0.2-py3-none-any.whl
# Creating
WORKDIR /app
COPY ./cellfree/prepare/protocol_cyclomics/tag_bam_by_json_metadata.py /app/
