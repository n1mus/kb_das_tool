/*
A KBase module: kb_das_tool
*/

module kb_das_tool {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef run_kb_das_tool(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;

};
