FROM python:3.11
RUN mkdir -p /signals/src
WORKDIR /signals
COPY setup.py setup.py
RUN python -m pip install -e .
COPY . .
ENTRYPOINT ["python", "-m", "signals"]
