# prop print

    Code
      print(prop(nFc > 0))
    Message
      prop(nFc > 0, na_action = "keep")

---

    Code
      print(prop(nFc > 0, within = T == "complex"))
    Message
      prop(nFc > 0, within = (T == "complex"), na_action = "keep")

---

    Code
      print(prop(nFc > 0, within = (T == "complex")))
    Message
      prop(nFc > 0, within = (T == "complex"), na_action = "keep")

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
      print(ratio(T == "complex", T == "hybrid"))
    Message
      ratio(T == "complex", T == "hybrid", na_action = "keep")

---

    Code
      print(ratio(T == "complex", T == "hybrid", within = T == "complex"))
    Message
      ratio(T == "complex", T == "hybrid", within = (T == "complex"), na_action =
      "keep")

---

    Code
      print(ratio(T == "complex", T == "hybrid", within = (T == "complex")))
    Message
      ratio(T == "complex", T == "hybrid", within = (T == "complex"), na_action =
      "keep")

---

    Code
      print(ratio(T == "complex", T == "hybrid", within = NULL))
    Message
      ratio(T == "complex", T == "hybrid", na_action = "keep")

---

    Code
      print(ratio(T == "complex", T == "hybrid", na_action = "zero"))
    Message
      ratio(T == "complex", T == "hybrid", na_action = "zero")

# wmean print

    Code
      print(wmean(nA))
    Message
      wmean(nA, na_action = "keep")

---

    Code
      print(wmean(nA, within = T == "complex"))
    Message
      wmean(nA, within = (T == "complex"), na_action = "keep")

---

    Code
      print(wmean(nA, within = (T == "complex")))
    Message
      wmean(nA, within = (T == "complex"), na_action = "keep")

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

