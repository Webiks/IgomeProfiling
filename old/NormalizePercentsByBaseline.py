
def normalize(percents): #with a name that is related to the file's name
    import numpy as np
    with open() as f:
        ncols = len(f.readline().split(','))

    data = np.loadtxt(percents, delimiter=',', skiprows=1, usecols=range(1, ncols + 1))

    #WAS DONE MANUALLY IN EXCEL
    #averaging the 3 baseline samples and dividing each sample by the baseline.


if __name__ == '__main__':
    import logging
    logger = logging.getLogger('main') # use logger instead of printing
    from sys import argv
    # This block will be executed only when you run it as your main program.
    # If this module is being imported from another script, this block won't be executed, however the function will be available...
    if len(argv) < 4: # change the number of arguments according to your script
        logger.error('Usage: python ' + argv[0] + ' table')
        exit()
    else:
        normalize(*argv[1:])
