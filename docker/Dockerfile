#
# Debian-based build and runtime enviroment for GreenALM
#
FROM gcc:9

ARG DEBIAN_FRONTEND=noninteractive

ARG USERNAME=developer
ARG USER_UID=1000
ARG USER_GID=$USER_UID

RUN echo "Acquire::http::Pipeline-Depth 0;" >> /etc/apt/apt.conf.d/99fixbadproxy && \
    echo "Acquire::http::No-Cache true;" >> /etc/apt/apt.conf.d/99fixbadproxy && \
    echo "Acquire::BrokenProxy    true;" >> /etc/apt/apt.conf.d/99fixbadproxy && \
    apt-get update && \
    apt-get install -y --no-install-recommends apt-utils && \
    apt-get upgrade -y && \
    apt-get install -y \
            curl \
            git \
            libarchive-dev \
            libffi-dev \
            libjsoncpp-dev \
            libgdbm-dev \
            libfftw3-dev \
            libhdf5-dev \
            liblapack-dev \
            liblzma-dev \
            libncurses-dev \
            libreadline-dev \
            librhash-dev \
            libsqlite3-dev \
            libssl-dev \
            libz-dev \
            ninja-build \
            procps \
            tk-dev \
            zlib1g && \
        apt-get autoremove -y && \
        apt-get clean -y && \
        rm -rf /var/lib/apt/lists/* && \
        mkdir -p /opt/build && \
   groupadd --gid $USER_GID $USERNAME && \
   useradd --uid $USER_UID --gid $USER_GID -m $USERNAME && \
   apt-get update && \
   apt-get install -y sudo && \
   echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME && \
   chmod 0440 /etc/sudoers.d/$USERNAME


WORKDIR /opt/build


#
# Install: OpenMPI
#
ARG OPENMPI_VERSION="4.1.4"
ARG OPENMPI_MAJOR_VERSION="v4.1"
ARG OPENMPI_CONFIGURE_OPTIONS
ARG OPENMPI_MAKE_OPTIONS="-j"

ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

RUN wget https://download.open-mpi.org/release/open-mpi/${OPENMPI_MAJOR_VERSION}/openmpi-${OPENMPI_VERSION}.tar.gz && \
    tar xfz openmpi-${OPENMPI_VERSION}.tar.gz && \
    cd openmpi-${OPENMPI_VERSION} && \
    ./configure ${OPENMPI_CONFIGURE_OPTIONS} && \
    make all ${OPENMPI_MAKE_OPTIONS} && \
    make install && \
    ldconfig && \
    mpicc --version

#
# Install: CMake
#
ARG CMAKE_VERSION="3.25.3"
ARG CMAKE_MAKE_OPTIONS="-j"

RUN curl -OL https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-SHA-256.txt && \
    curl -OL https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}.tar.gz && \
    tar -xf cmake-${CMAKE_VERSION}.tar.gz && \
    cd cmake-${CMAKE_VERSION} && \
    ./bootstrap && \
    make ${CMAKE_MAKE_OPTIONS} && \
    make install && \
    cmake --version


#
# Install: LibXC
#
ARG LIBXC_VERSION="5.2.3"
ARG LIBXC_MAKE_OPTIONS="-j"

RUN git clone -b ${LIBXC_VERSION} --depth 1 --quiet  https://gitlab.com/libxc/libxc.git && \
    cd libxc && \
    mkdir -p build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=Debug \
          -DNAMESPACE_INSTALL_INCLUDEDIR=/  \
          -DBUILD_SHARED_LIBS=ON \
          -DENABLE_FORTRAN=ON \
          -DBUILD_FPIC=ON \
          -DBUILD_TESTING=ON .. && \
    make ${LIBXC_MAKE_OPTIONS} && make test && make install && \
    ldconfig

#
# Zsh
#
RUN sh -c "$(wget -O- https://github.com/deluan/zsh-in-docker/releases/download/v1.1.2/zsh-in-docker.sh)" -- \
    -p git \
    -p https://github.com/zsh-users/zsh-autosuggestions \
    -p https://github.com/zsh-users/zsh-completions \
    -p https://github.com/zsh-users/zsh-history-substring-search \
    -p https://github.com/zsh-users/zsh-syntax-highlighting \
    -p 'history-substring-search' \
    -a 'bindkey "\$terminfo[kcuu1]" history-substring-search-up' \
    -a 'bindkey "\$terminfo[kcud1]" history-substring-search-down'

ENV LIBXC_LIBRARIES=/usr/local/lib/libxc.so
ENV LIBXC_INCLUDE_DIR=/usr/local/include/



# User configurations
ENV PYTHONUNBUFFERED=1

ENTRYPOINT [ "/bin/zsh" ]
CMD [ "-l" ]



