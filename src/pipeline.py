import argparse
import sys
import parameters as param

import hybrid
import polish
#import graph

COMMAND = "command"
HYBRID  = "hybrid"
POLISH  = "polish"
GRAPH   = "graph"

if __name__== "__main__":
    
    mainParser = argparse.ArgumentParser(description="Hybrid assembly tool")
    subparsers = mainParser.add_subparsers(dest=COMMAND)
    
    hybridParser = subparsers.add_parser(HYBRID)
    polishParser = subparsers.add_parser(POLISH)
    graphParser = subparsers.add_parser(GRAPH)
    
    hybridParser = param.set_hybrid_parameters(hybridParser)
    polishParser = param.set_polish_parameters(polishParser)
    graphParser = param.set_graph_parameters(graphParser)

    command = vars(mainParser.parse_args())[COMMAND]
    if command is not None:

        if command == HYBRID:
            #parameters = param.get_hybrid_parameters(hybridParser)
            parameters = param.explicit_hybrid_params()

            if not param.validate_hybrid_parameters(parameters):
                hybridParser.print_help()
            else:
                hybrid.main(parameters)
            sys.exit()

        if command == POLISH:
            #parameters = param.get_polish_parameters(polishParser)
            parameters = param.explicit_polish_params()

            if not param.validate_polish_parameters(parameters):
                polishParser.print_help()
            else:
                polish.main(parameters)
            sys.exit()

        if command == GRAPH:
            parameters = param.get_polish_parameters(graphParser)
            if not param.validate_graph_parameters(parameters):
                graphParser.print_help()
            else:
                #graph.main(parameters)
                pass
            sys.exit()

    mainParser.print_help()
    sys.exit()
