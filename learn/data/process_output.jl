function get_performance(filename)
    file = open(filename)
    old_line = ""
    accs = Float64[]
    for line in eachline(file)
      if (line[1:3] == "Hit")
        temp = parse(split(old_line)[6])
        push!(accs, temp)
      end
      old_line = line
    end
    return accs
end


function main()
  num_args = 1
  filename = ARGS[num_args]; num_args += 1;
  c = get_performance(filename)
  println(c) 
  println(mean(c))
  println(std(c))
end


main()
