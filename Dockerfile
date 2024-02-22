FROM ubuntu:22.04

# Install system dependencies
RUN apt-get update && apt-get install -y build-essential python3-pip virtualenv \
    python3-venv sqlite3 libsqlite3-dev sqlcipher libsqlcipher-dev git

# Set user, group, and create home directory
ARG user=appuser
ARG group=appuser
ARG uid=1000
ARG gid=1000
RUN groupadd -g ${gid} ${group}
RUN useradd -u ${uid} -g ${group} -s /bin/sh -m ${user}

# Aliases for application and working directories
ENV appdir=/home/${user}/app
ENV workdir=/home/${user}/work

WORKDIR ${appdir}

# Copy project files to the application directory
COPY aware aware
COPY pyproject.toml .
COPY MANIFEST.in .

ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

# Create a new virtual environment, install dependencies, and AWARE
RUN python3 -m venv .venv && . .venv/bin/activate && python -m pip install -U pip wheel setuptools .

# Create entrypoint for application
WORKDIR ${workdir}

ENTRYPOINT [ "/bin/bash", "-c", "exec ${appdir}/.venv/bin/python -m aware \"${@}\"", "--" ]

# Optional arguments (-t means run Telegram bot)

CMD ["-t"]
