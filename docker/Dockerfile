FROM python:3
#ADD index.html index.html
RUN python -m pip install --upgrade pip
RUN python -m pip install fastapi uvicorn brdr
ADD grb_webservice.py grb_webservice.py
EXPOSE 7999
ENTRYPOINT ["python3", "grb_webservice.py"]