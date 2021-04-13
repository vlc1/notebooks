using Pluto

function pluto2html(filename)
	ss = Pluto.ServerSession()
	nb = Pluto.SessionActions.open(ss, filename; run_async=false)
	html = Pluto.generate_html(nb)
	write(replace(filename, r".jl"i => ".html"), html)
end

