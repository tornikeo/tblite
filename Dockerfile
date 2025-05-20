# Use an official NVIDIA CUDA base image with Python pre-installed
FROM nvidia/cuda:11.8.0-devel-ubuntu22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

# Install required dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    liblapack-dev \
    build-essential \
    gfortran \
    python3 \
    python3-pip \
    python3-dev \
    ninja-build \
    cmake \
    git \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install PyTorch
RUN pip3 install --no-cache-dir torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

# Set working directory
WORKDIR /workspace

# Copy the project files into the container
COPY . /workspace

RUN pip install --no-cache-dir meson ninja ipykernel pandas matplotlib scikit-learn

# Build the project using Meson
RUN pip3 install meson && \
    meson setup _build && \
    meson compile -C _build

# Default command
CMD ["bash"]