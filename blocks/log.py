def out(text, verboseLevel, param, wait=False):
    if param.VERBOSE >= verboseLevel:
        print(text)
        if wait:
            input()
    
