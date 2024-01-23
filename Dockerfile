########################################################################################################################
# build stage
########################################################################################################################


FROM alpine:3.15.10 AS build
RUN apk update && \
    apk add --no-cache \
        build-base \
        cmake \
        git \
        gcc \
        g++ \
        libstdc++ \
        libgomp

WORKDIR /taxor

COPY src/ ./src/

WORKDIR /taxor/build

RUN cmake ../src && \
    cmake --build . --config Release

########################################################################################################################
# image
########################################################################################################################

FROM alpine:3.15.10

RUN apk update && \
    apk add --no-cache \
    libstdc++ \
    libgomp

RUN addgroup -S taxor && adduser -S taxor -G taxor
USER taxor

COPY --chown=taxor:taxor --from=build \
    ./taxor/build/main/taxor \
    ./app/

ENTRYPOINT [ "./app/taxor" ]