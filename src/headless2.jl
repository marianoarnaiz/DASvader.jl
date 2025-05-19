"A file that is a copy of headless but not part of the module"

function env_has_display()
    if Sys.islinux()
        return !isempty(get(ENV, "DISPLAY", "")) || !isempty(get(ENV, "WAYLAND_DISPLAY", ""))
    elseif Sys.iswindows() || Sys.isapple()
        return true  # Assume GUI unless proven otherwise
    else
        return false
    end
end

function glmakie_screen_available()
    try
        @eval using GLMakie  # Delay loading to avoid errors
        s = GLMakie.Screen()
        close(s)
        return true
    catch
        return false
    end
end

function has_graphical_display()
    env_has_display() && glmakie_screen_available()
end

export has_graphical_display, glmakie_screen_available, env_has_display
