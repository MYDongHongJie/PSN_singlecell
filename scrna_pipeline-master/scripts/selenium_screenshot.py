#!/usr/bin/env python3
# encoding: utf-8
# Author: luyao
# Desc: This script is used to take screenshots automatically.

import os, time, argparse
from selenium import webdriver
from PIL import Image

#=======================================================================================================================
parser = argparse.ArgumentParser(description="质控自动截图")
parser.add_argument("-i", "--input", help="web_summary.html")
parser.add_argument("-n", "--name", help="the name of output file")
parser.add_argument("-o", "--outdir", default="result/cellranger", help="the output directory of screenshot")
args = parser.parse_args()
html = os.path.abspath(args.input)
filename = args.name
outdir = os.path.abspath(args.outdir)

def screenshot(html, filename, outdir):
    '''
    :param html:
    :param outdir:
    :return: screenshot of web_summary.html
    '''
    options = webdriver.ChromeOptions()
    options.add_argument("headless")
    options.add_argument('no-sandbox')
    options.add_argument('dns-prefetch-disable')
    options.add_argument('disable-gpu')
    options.add_argument('disable-dev-shm-usage')
    options.add_argument('disable-features=VizDisplayCompositor')
    options.add_argument('disable-features=NetworkService')
    options.add_argument('window-size=1920x945')
    options.add_experimental_option("excludeSwitches",["ignore-certificate-errors"])
    driver = webdriver.Chrome('chromedriver', chrome_options=options)
    driver.set_page_load_timeout(10)
    driver.get('file://' + html)
    driver.set_window_size(1920, 945)
    driver.maximize_window()  # full screen
    width = driver.execute_script("return document.documentElement.scrollWidth")
    height = driver.execute_script("return document.documentElement.scrollHeight")
    # print(width,height)
    driver.set_window_size(width, height)
    # driver.get_window_size()  # Get the size of the current window
    # time.sleep(1)
    driver.save_screenshot(outdir + "/" + filename + ".png")
    driver.close()
    # Resize the picture
    img = Image.open(outdir + "/" + filename + ".png")
    width, height = img.size
    box = (300, 0, 1550, height)
    region = img.crop(box)
    region.save(outdir + "/" + filename + ".png")

if __name__ == '__main__':
    screenshot(html, filename, outdir)
