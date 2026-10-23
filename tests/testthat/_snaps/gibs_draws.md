# process_I_multivariant()

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["local", "imported"]
        }
      },
      "value": [
        {
          "type": "double",
          "attributes": {},
          "value": [10, 10, 10, 10, 10]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [0, 0, 0, 0, 0]
        }
      ]
    }

# compute_lambda()

    {
      "type": "double",
      "attributes": {},
      "value": [10, 10, 10, 10, 10]
    }

# draw_R()

    {
      "type": "double",
      "attributes": {},
      "value": [1.34883043, 0.84282641, 0.79784422, 1.38220267, 0.78447105]
    }

# estimate_advantage()

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
          "value": [1.00471817, 1.03107137, 1.02035086, 1.00506093, 1.00732389]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [0.91800364, 0.85016542, 0.71860597, 1.07449868, 1.30864555]
        },
        {
          "type": "logical",
          "attributes": {},
          "value": [false]
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
                  "value": [1.02097152, 1.11993594]
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

