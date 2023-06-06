FROM python:3.11
RUN mkdir -p /nodclust/src
WORKDIR /nodclust
COPY setup.py setup.py
RUN python -m pip install -e .
COPY . .
ENTRYPOINT ["python", "-m", "nodclust"]
