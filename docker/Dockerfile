FROM python:3
#ADD index.html index.html
RUN python -m pip install --upgrade pip
RUN python -m pip install pyramid brdr shapely
ADD grb_webservice.py grb_webservice.py
EXPOSE 8000
ENTRYPOINT ["python3", "grb_webservice.py"]