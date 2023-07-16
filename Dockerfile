FROM python:3.11
RUN mkdir -p /modena/src
WORKDIR /modena
COPY setup.py setup.py
RUN python -m pip install --upgrade pip
RUN python -m pip install -e .
COPY . .
ENTRYPOINT ["python", "-m", "modena"]
