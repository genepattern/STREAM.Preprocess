FROM pinellolab/stream:0.3.8

RUN mkdir /stream
COPY preprocess_command_line.py /stream/preprocess_command_line.py

ENTRYPOINT []
