#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: Score predictions file

requirements:
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
    - entryname: score.py
      entry: |
        #!/usr/bin/env python
        import argparse
        import json
        parser = argparse.ArgumentParser()
        parser.add_argument("-f", "--submissionfile", required=True, help="Submission File")
        parser.add_argument("-r", "--results", required=True, help="Scoring results")
        parser.add_argument("-g", "--goldstandard", required=True, help="Goldstandard for scoring")

        args = parser.parse_args()

        import yaml
        import numpy as np
        from scoring import ScoreSweeps

        sweepscore = ScoreSweeps(args.submissionfile, args.goldstandard)

        sweepscore.calculate_stats()
        true_pos_prop = sweepscore.true_positive_score
        false_neg_prop = sweepscore.false_negative_score
        f1_score = sweepscore.f1

        prediction_file_status = "SCORED"

        result = {'true_positive_score': true_pos_prop,
                  'false_negative_score': false_neg_prop,
                  'f1_score': f1_score,
                  'submission_status': prediction_file_status}
        with open(args.results, 'w') as o:
          o.write(json.dumps(result))

inputs:
  - id: input_file
    type: File
  - id: goldstandard
    type: File
  - id: check_validation_finished
    type: boolean?

outputs:
  - id: results
    type: File
    outputBinding:
      glob: results.json
  - id: status
    type: string
    outputBinding:
      glob: results.json
      outputEval: $(JSON.parse(self[0].contents)['submission_status'])
      loadContents: true

baseCommand: python
arguments:
  - valueFrom: score.py
  - prefix: -f
    valueFrom: $(inputs.input_file.path)
  - prefix: -g
    valueFrom: $(inputs.goldstandard.path)
  - prefix: -r
    valueFrom: results.json

hints:
  DockerRequirement:
    dockerPull: tjstruck/popsim-pilot-slim:1.32
