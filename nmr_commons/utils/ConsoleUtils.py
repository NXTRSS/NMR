import sys
import numpy as np

# Print iterations progress
def printProgress (iteration, total, prefix = '', suffix = '', decimals = 2, barLength = 100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : number of decimals in percent complete (Int)
        barLength   - Optional  : character length of bar (Int)
    """
    filledLength    = int(round(barLength * iteration / float(total)))
    percents        = round(100.00 * (iteration / float(total)), decimals)
    bar             = '|' * filledLength + '-' * (barLength - filledLength)
    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),
    sys.stdout.flush()
    if iteration == total:
        sys.stdout.write('\n')
        sys.stdout.flush()


def to_nice_string(iterable, form=None, brackets=None, sep=' '):
    """
    Args:
        iterable: list, tuple, set, ndarray or any iterable collection
        form: format string for single element, e.g. '2.3f' or '<5'
        brackets: brackets to use at the beginning and at the end of the string, e.g. '|' or '()'
        sep: separator

    Returns:
        nicely formatted string
    """
    # TODO: handle form=None (intelligent format selection)
    if form is None:
        pass

    form = '{{:{}}}'.format(form)
    string = sep.join([form.format(elem) for elem in iterable])
    if brackets is not None:
        if len(brackets) == 2:
            string = brackets[0] + string + brackets[1]
        else:
            string = brackets + string + brackets
    return string
