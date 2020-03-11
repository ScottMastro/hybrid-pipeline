class fail_enum:
    #failure enums
    FULL_OVERLAP      = 1
    PARTIAL_OVERLAP   = 2
    EXCESSIVE_TRIM    = 3
    SAME_PATH         = 4
    PEEL_FAIL         = 5
    NO_UNITIG         = 6
    PASS              = 0
    UNSURE            = -1
    
    failDict = {UNSURE:"????", PASS:"pass", FULL_OVERLAP:"overlap (full)",
                PARTIAL_OVERLAP:"overlap (partial)", 
                EXCESSIVE_TRIM:"excessive trimming", SAME_PATH:"redundant",
                NO_UNITIG:"no unitig support", PEEL_FAIL:"could not find critical forks"}
