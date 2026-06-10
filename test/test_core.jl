using FiniteDiff
using ForwardDiff
using Dates
using TimeZones

@testset "PVlib basic pipeline" begin
    @testset "AC power end-to-end (deterministic)" begin
        # Representative midday weather sample for a single-site regression.
        weather_data = WeatherSample{Float64}(
            ZonedDateTime(2020, 6, 1, 12, 0, 0, TimeZone("America/Denver")),
            800.0,
            900.0,
            100.0,
            20.0,
            20.0,
            1.0,
            101325.0,
            1.0,
            1.0,
            0.1,
        )

        # Site location and fixed-tilt array geometry.
        latitude = 35.1
        longitude = -106.6
        altitude = 1500.0
        surface_tilt = latitude
        surface_azimuth = 180.0
        albedo = 0.25

        # Run the explicit PV pipeline from solar position through AC power.
        solar_pos = get_solar_position(latitude, longitude, altitude, weather_data)
        total_irradiance = get_total_irradiance(
            surface_tilt,
            surface_azimuth,
            weather_data,
            solar_pos,
            albedo,
        )

        pv_module = read_solar_module("Canadian Solar CS5P-220M [ 2009]")
        pv_inverter = read_solar_inverter("ABB: MICRO-0.25-I-OUTD-US-208 [208V]")

        cell_temp = sapm_cell_temperature(total_irradiance, weather_data)
        effective_irradiance = sapm_effective_irradiance(
            total_irradiance,
            pv_module,
            solar_pos,
            surface_tilt,
            surface_azimuth,
            altitude,
        )

        dc_power = sapm_dc_components(pv_module, effective_irradiance, cell_temp)
        ac_power = sandia_ac_power(pv_inverter, dc_power)

        @test isapprox(ac_power.ac_power, 167.58832880266323; atol = 1e-6)

        # Cover vector overloads for tilt-only, azimuth-only, and combined dispatch.
        weather_vec = [weather_data, weather_data]
        solar_vec = [solar_pos, solar_pos]
        tilt_vec = [surface_tilt, surface_tilt + 5.0]
        azimuth_vec = [surface_azimuth, surface_azimuth - 10.0]

        irr_tilt =
            get_total_irradiance(tilt_vec, surface_azimuth, weather_vec, solar_vec, albedo)
        irr_az =
            get_total_irradiance(surface_tilt, azimuth_vec, weather_vec, solar_vec, albedo)
        irr_both =
            get_total_irradiance(tilt_vec, azimuth_vec, weather_vec, solar_vec, albedo)

        @test length(irr_tilt) == 2
        @test length(irr_az) == 2
        @test length(irr_both) == 2
        @test all(isfinite(getfield(x, :poa_global)) for x in irr_both)

        # Compare AD and finite differencing through the full AC-power path.
        tilt_to_ac(tilt) = begin
            irr = get_total_irradiance(tilt, surface_azimuth, weather_data, solar_pos, albedo)
            cell_temp_local = sapm_cell_temperature(irr, weather_data)
            eff = sapm_effective_irradiance(
                irr,
                pv_module,
                solar_pos,
                tilt,
                surface_azimuth,
                altitude,
            )
            dc = sapm_dc_components(pv_module, eff, cell_temp_local)
            sandia_ac_power(pv_inverter, dc).ac_power
        end

        dP_dtilt_ad = ForwardDiff.derivative(tilt_to_ac, surface_tilt)
        dP_dtilt_fd =
            FiniteDiff.finite_difference_derivative(tilt_to_ac, surface_tilt, Val(:central))
        @test isfinite(dP_dtilt_ad)
        @test isfinite(dP_dtilt_fd)
        @test isapprox(dP_dtilt_ad, dP_dtilt_fd; rtol = 1.0e-5, atol = 1.0e-7)

        # Compare AD and finite differencing for the irradiance-only sensitivity.
        ghi_to_poa(ghi) = begin
            weather = WeatherSample{typeof(ghi)}(
                weather_data.time,
                ghi,
                weather_data.dni,
                weather_data.dhi,
                weather_data.temp_air,
                weather_data.temp_dew,
                weather_data.relative_humidity,
                weather_data.pressure,
                weather_data.wind_speed,
                weather_data.wind_direction,
                weather_data.albedo,
            )
            get_total_irradiance(surface_tilt, surface_azimuth, weather, solar_pos, albedo).poa_global
        end

        dpoa_dghi_ad = ForwardDiff.derivative(ghi_to_poa, weather_data.ghi)
        dpoa_dghi_fd = FiniteDiff.finite_difference_derivative(
            ghi_to_poa,
            weather_data.ghi,
            Val(:central),
        )
        @test isfinite(dpoa_dghi_ad)
        @test isfinite(dpoa_dghi_fd)
        @test dpoa_dghi_ad > 0
        @test isapprox(dpoa_dghi_ad, dpoa_dghi_fd; rtol = 1.0e-6, atol = 1.0e-8)
    end
end
