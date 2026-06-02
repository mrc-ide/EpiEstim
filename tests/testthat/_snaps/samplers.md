# estimate_advantage produces expected results (2 var, 2 loc, R_loc1 = 1.1, R_loc2 = 1.5)

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["epsilon", "R", "convergence", "diag"]
        }
      },
      "value": [
        {
          "type": "double",
          "attributes": {},
          "value": [1.500036, 1.500023, 1.500019, 1.50001, 1.500015]
        },
        {
          "type": "double",
          "attributes": {},
          "value": ["NA", 1.500005, 1.013014, 1.499871, 1.10004]
        },
        {
          "type": "logical",
          "attributes": {},
          "value": [true]
        },
        {
          "type": "list",
          "attributes": {},
          "value": [
            {
              "type": "list",
              "attributes": {
                "names": {
                  "type": "character",
                  "attributes": {},
                  "value": ["psrf", "mpsrf"]
                }
              },
              "value": [
                {
                  "type": "double",
                  "attributes": {
                    "dim": {
                      "type": "integer",
                      "attributes": {},
                      "value": [1, 2]
                    },
                    "dimnames": {
                      "type": "list",
                      "attributes": {},
                      "value": [
                        {
                          "type": "NULL"
                        },
                        {
                          "type": "character",
                          "attributes": {},
                          "value": ["Point est.", "Upper C.I."]
                        }
                      ]
                    }
                  },
                  "value": [0.997933, 1.056718]
                },
                {
                  "type": "NULL"
                }
              ]
            }
          ]
        }
      ]
    }

