# t_test

    Code
      t_test(efc, "e17age")
    Output
      # One Sample t-test
      
        Data: e17age
        Group 1: e17age (mean = 79.12)
        Alternative hypothesis: true mean is not equal to 0
      
        t = 291.78, Cohen's d = 9.77 (large effect), df = 890, p < .001
      

---

    Code
      t_test(efc, "e17age", "e16sex")
    Output
      # Welch Two Sample t-test
      
        Data: e17age by e16sex
        Group 1: 1 (n = 294, mean = 76.16)
        Group 2: 2 (n = 596, mean = 80.57)
        Alternative hypothesis: true difference in means is not equal to 0
      
        t = -8.05, Cohen's d = -0.56 (medium effect), df = 610.8, p < .001
      

---

    Code
      t_test(efc, c("e17age", "c160age"))
    Output
      # Welch Two Sample t-test
      
        Data: e17age by c160age
        Group 1: c160age (n = 890, mean = 53.42)
        Group 2: e17age  (n = 890, mean = 79.12)
        Alternative hypothesis: true difference in means is not equal to 0
      
        t = -49.22, Cohen's d = -2.33 (large effect), df = 1468.1, p < .001
      

---

    Code
      t_test(efc, c("e17age", "c160age"), paired = TRUE)
    Output
      # Paired t-test
      
        Data: e17age and c160age (mean difference = 25.70)
        Alternative hypothesis: true mean is not equal to 0
      
        t = 54.11, Cohen's d = 1.81 (large effect), df = 889, p < .001
      

---

    Code
      t_test(efc, "e17age", weights = "weight")
    Output
      # One Sample t-test (weighted)
      
        Data: e17age
        Group 1: e17age (n = 897, mean = 79.17)
        Alternative hypothesis: true mean is not equal to 0
      
        t = 291.31, Cohen's d = 3.17 (large effect), df = 890, p < .001
      

---

    Code
      t_test(efc, "e17age", "e16sex", weights = "weight")
    Output
      # Two-Sample t-test (weighted)
      
        Data: e17age by e16sex
        Group 1: 1 (n = 600, mean = 80.63)
        Group 2: 2 (n = 296, mean = 76.19)
        Alternative hypothesis: true difference in means is not equal to 0
      
        t = 8.03, Cohen's d = -0.17 (very small effect), df = 604.5, p < .001
      

---

    Code
      t_test(efc, c("e17age", "c160age"), weights = "weight")
    Output
      # Two-Sample t-test (weighted)
      
        Data: e17age by c160age
        Group 1: c160age (n = 896, mean = 79.17)
        Group 2: e17age  (n = 896, mean = 53.40)
        Alternative hypothesis: true difference in means is not equal to 0
      
        t = 49.31, Cohen's d = -1.12 (large effect), df = 1470.0, p < .001
      

---

    Code
      t_test(efc, c("e17age", "c160age"), weights = "weight", paired = TRUE)
    Output
      # Paired t-test (weighted)
      
        Data: e17age and c160age (mean difference = 25.77)
        Alternative hypothesis: true mean difference is not equal to 0
      
        t = 54.37, Cohen's d = 1.54 (large effect), df = 889, p < .001
      

