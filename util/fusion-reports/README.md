

## Creating a fusion inspector report 

Note: These instructions are specifically for Fusion Inspector reports.  See [https://github.com/FusionInspector/FusionInspector/wiki](https://github.com/FusionInspector/FusionInspector/wiki)

Fusion inspector reports are created by combining Fusion Inpsector output with an html template file.

First make sure that the html template file contains the comment line `<!-- start igv report here -->` within a script tag.
Directly below this line a javscript variable called "data" will be created so insure that no other variable names conflict.  
  
Then use the create_fusion_report.py script to create a new self-contained html file.
```sh
python create_fusion_report.py filename fusions
```
where filename is the path to the igv.js html template file and fusions is a fusion inspector output file defining
the fusions

To run the example execute

```sh
python create_fusion_report.py examples/fusions/igvjs_fusion.html
```

After, running the script, see examples/fusions/igvjs_fusion_viewer_report.html for the result.


