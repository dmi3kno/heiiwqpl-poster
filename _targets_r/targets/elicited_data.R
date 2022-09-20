tar_target(elicited_data, {
  tibble::tribble(~cat_idx, ~cat, ~CDFprob, ~count,
                    1L, "Morning coffee (0-4 g/day)", 0.25, 20,#18
                    1L, "Morning coffee (0-4 g/day)", 0.5, 24,#24
                    1L, "Morning coffee (0-4 g/day)", 0.75, 28,#30
                    2L, "After dinner (4-8 g/day)", 0.25, 49-24, #49
                    2L, "After dinner (4-8 g/day)", 0.50, 55-24, #55
                    2L, "After dinner (4-8 g/day)", 0.75, 60-24, #60
                    4L, "Sweet tooth (16+ g/day)", 0.25, 100-95,#92
                    4L, "Sweet tooth (16+ g/day)", 0.50, 100-85,#85
                    4L, "Sweet tooth (16+ g/day)", 0.75, 100-80,#80
    )
})
