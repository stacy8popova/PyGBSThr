FROM archlinux:latest
RUN pacman -Sy --noconfirm --noprogressbar \
        gcc \
        make \
        python \
        python-pip \
        unzip

RUN pip install numpy matplotlib sympy thewalrus strawberryfields jupyterlab
ADD ./ /home/devel/

USER root
WORKDIR /home/devel
ENV HOME=/home/devel
RUN make
CMD ["jupyter", "notebook", "--port=8080", "--no-browser", "--ip=0.0.0.0", "--allow-root"]


#docker run -p 8080:8080 stacy8popova/py_gbs_thr

#docker run -it -v $(pwd):/home/devel/in_out stacy8popova/py_gbs_thr python3 ./start_up.py
#docker run -it -v ${pwd}:/home/devel/in_out stacy8popova/py_gbs_thr python3 ./start_up.py
#docker run -it -v %cd%:/home/devel/in_out stacy8popova/py_gbs_thr python3 ./start_up.py