function matrix2table(filename::String, outputname::String)

  # Spacer
  s = repeat(" ", 2)

  # Load matrix and format as string
  body = Vector{String}()
  open(filename) do file
    for ln in eachline(file)
      push!(body, join([replace(ln, "\t" => " & "), "\\\\"]))
    end
  end
  header = body[1]
  body = body[2:end]
  num_cols = count(i -> (i == '&'), body[1]) + 1

  # Create output .tex file
  io = open(outputname, "w")
    #
    write(io, "% This file was generated by Julia2pgf.\n")
    write(io, "\\begin{tabularx}{\\textwidth}{$(repeat("X",num_cols))}\n")
    write(io, "$(s)\\toprule\n")
    write(io, "$(s)$(header)\n")
    write(io, "$(s)\\midrule\n")
    for ln in body
      write(io, "$(s)$(ln)\n")
    end
    write(io, "$(s)\\bottomrule\n")
    write(io, "\\end{tabularx}\n")
  close(io)
end
