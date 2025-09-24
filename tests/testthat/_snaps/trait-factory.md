# prop print

    Code
      print(prop(nFc > 0))
    Message
      prop(nFc > 0, na_action = "keep")

---

    Code
      print(prop(nFc > 0, within = Tp == "complex"))
    Message
      prop(nFc > 0, within = (Tp == "complex"), na_action = "keep")

---

    Code
      print(prop(nFc > 0, within = (Tp == "complex")))
    Message
      prop(nFc > 0, within = (Tp == "complex"), na_action = "keep")

---

    Code
      print(prop(nFc > 0, within = NULL))
    Message
      prop(nFc > 0, na_action = "keep")

---

    Code
      print(prop(nFc > 0, na_action = "zero"))
    Message
      prop(nFc > 0, na_action = "zero")

# ratio print

    Code
      print(ratio(Tp == "complex", Tp == "hybrid"))
    Message
      ratio(Tp == "complex", Tp == "hybrid", na_action = "keep")

---

    Code
      print(ratio(Tp == "complex", Tp == "hybrid", within = Tp == "complex"))
    Message
      ratio(Tp == "complex", Tp == "hybrid", within = (Tp == "complex"), na_action =
      "keep")

---

    Code
      print(ratio(Tp == "complex", Tp == "hybrid", within = (Tp == "complex")))
    Message
      ratio(Tp == "complex", Tp == "hybrid", within = (Tp == "complex"), na_action =
      "keep")

---

    Code
      print(ratio(Tp == "complex", Tp == "hybrid", within = NULL))
    Message
      ratio(Tp == "complex", Tp == "hybrid", na_action = "keep")

---

    Code
      print(ratio(Tp == "complex", Tp == "hybrid", na_action = "zero"))
    Message
      ratio(Tp == "complex", Tp == "hybrid", na_action = "zero")

# wmean print

    Code
      print(wmean(nA))
    Message
      wmean(nA, na_action = "keep")

---

    Code
      print(wmean(nA, within = Tp == "complex"))
    Message
      wmean(nA, within = (Tp == "complex"), na_action = "keep")

---

    Code
      print(wmean(nA, within = (Tp == "complex")))
    Message
      wmean(nA, within = (Tp == "complex"), na_action = "keep")

---

    Code
      print(wmean(nA, within = NULL))
    Message
      wmean(nA, na_action = "keep")

---

    Code
      print(wmean(nA, na_action = "zero"))
    Message
      wmean(nA, na_action = "zero")

