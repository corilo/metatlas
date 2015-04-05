"""Selected functions from https://github.com/kbase/transform/blob/master/lib/biokbase/Transform/script_utils.py
For use in Travis testing
"""

import logging
import re
import sys
import time


def stderrlogger(name, level=logging.INFO):
    """
    Return a standard python logger with a stderr handler attached and using a prefix
    format that will make logging consistent between scripts.
    """

    logger = logging.getLogger(name)
    logger.setLevel(level)

    # send messages to sys.stderr
    streamHandler = logging.StreamHandler(sys.stderr)

    formatter = logging.Formatter("%(asctime)s - %(filename)s - %(lineno)d - %(levelname)s - %(message)s")
    formatter.converter = time.gmtime
    streamHandler.setFormatter(formatter)

    logger.addHandler(streamHandler)

    return logger


def parse_docs(docstring=None):
    """
    Parses the docstring of a function and returns a dictionary of the elements.
    """

    # TODO, revisit this, probably can use other ways of doing this
    script_details = dict()

    keys = ["Authors","Returns","Args"]

    remainder = docstring[:]
    for k in keys:
        remainder, script_details[k] = remainder.split(k+":",1)
        script_details[k] = script_details[k].strip()

    script_details["Description"] = remainder

    # special treatment for Args since we want a dict, split on :, then cleanup whitespace
    # keep the : in the keys to do clean splits when getting the values
    argument_keys = [x.strip() for x in re.findall(".*:",script_details["Args"])]

    # split on the string in reverse by the keys, then wash out the extra whitespace
    remainder = script_details["Args"]
    argument_values = list()
    for k in reversed(argument_keys):
        remainder, value = remainder.split(k)
        argument_values.append(" ".join([x.strip() for x in value.split("\n")]))

    # create the dict using they keys without :, then get the values in the correct order
    script_details["Args"] = dict(zip([x.replace(":","") for x in argument_keys], reversed(argument_values)))

    return script_details



def upload_file_to_shock(logger = stderrlogger(__file__),
                         shock_service_url = None,
                         filePath = None,
                         ssl_verify = True,
                         token = None):
    """
    Use HTTP multi-part POST to save a file to a SHOCK instance.
    """

    if token is None:
        raise Exception("Authentication token required!")

    #build the header
    header = dict()
    header["Authorization"] = "Oauth {0}".format(token)

    if filePath is None:
        raise Exception("No file given for upload to SHOCK!")

    dataFile = open(os.path.abspath(filePath), 'r')
    m = MultipartEncoder(fields={'upload': (os.path.split(filePath)[-1], dataFile)})
    header['Content-Type'] = m.content_type

    logger.info("Sending {0} to {1}".format(filePath,shock_service_url))
