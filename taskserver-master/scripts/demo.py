import sys
import time

import fire


def demo(sleep: int, **kwargs):
    print('Hello OE Biotech!')
    print(sys.argv)
    print(kwargs)
    time.sleep(sleep)
    print('task finish~~~~~')


if __name__ == '__main__':
    fire.Fire()
