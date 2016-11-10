#! /bin/bash

for dm in 21.0 21.4 21.8 22.2 22.6 23.0 23.4 23.8 24.2 24.6 25.0 25.4 25.8 26.2 26.6 27.0
do
	python magfilter.py --iso --fwhm=2.0 --dm=$dm
done

PDFconcat --shuffle --output f_all_2.0.pdf f_*.pdf
