#version 330 core
out vec4 FragColor;

in vec4 sphere_obj_world; // Sphere center (xyz) and radius (w) in world space
in vec3 fpos_world;       // Fragment position on the bounding mesh in world space
in vec4 atomColor_out;    // Color from vertex shader

// Uniforms from BaseGLWidget
uniform vec3 viewPos;      // Camera position in world space
uniform mat4 projection;
uniform mat4 view;
uniform vec3 lightPos;
uniform vec3 lightColor;

// Lighting parameters (can be tuned or made into uniforms)
const vec3  ambient_material_color = vec3(0.2, 0.2, 0.2); // Base ambient color
const float specular_strength_factor = 0.7;
const float shininess_factor = 128.0;

vec2 rayPointDist( vec3 ray0, vec3 hRay, vec3 point ){
	vec3 pt  = point - ray0;	
	float t  = dot( hRay, pt );
    pt      -= t * hRay;
	float r2 = dot( pt, pt );  
	return vec2( t, r2 );
}

float raySphere( vec3 ray0, vec3 hRay, vec3 center, float R ){
	vec2 res = rayPointDist( ray0, hRay, center );
	float dr2 = R*R - res.y;
	if( dr2 > 0.0 ){
		return res.x - sqrt( dr2 );
	}else{
		return 1e+8; // No intersection or behind ray origin
	}
}

vec3 sphereNormal( float t, vec3 ray0, vec3 hRay, vec3 center ){
	return normalize(ray0 + t*hRay - center);
}

void main()
{
    vec3 ray_origin = viewPos; 
    vec3 ray_direction = normalize(fpos_world - viewPos); 

    vec3 sphere_center_w = sphere_obj_world.xyz;
    // sphere_radius_w will be the small radius (e.g., 0.05) from instanceActualSphereRadius
    float sphere_radius_w = sphere_obj_world.w; 

    float t = raySphere( ray_origin, ray_direction, sphere_center_w, sphere_radius_w );

    if( t > 1e+5 || t < 0.0 ){ // Discard if no hit, or hit behind camera
       discard;
    } else {
        vec3 P_intersect_world = ray_origin + t * ray_direction;
        vec3 N = sphereNormal( t, ray_origin, ray_direction, sphere_center_w );

        // Lighting
        vec3 ambient = ambient_material_color * lightColor;

        vec3 L = normalize(lightPos - P_intersect_world); // Light direction
        float diff_intensity = max(dot(N, L), 0.0);
        vec3 diffuse = diff_intensity * lightColor;

        vec3 V = normalize(viewPos - P_intersect_world); // View direction
        vec3 R = reflect(-L, N); // Reflected light direction
        float spec_intensity = pow(max(dot(V, R), 0.0), shininess_factor);
        vec3 specular = specular_strength_factor * spec_intensity * lightColor;

        vec3 result_rgb = (ambient + diffuse + specular) * atomColor_out.rgb;
        FragColor = vec4(result_rgb, atomColor_out.a);

        // Calculate gl_FragDepth
        vec4 clip_space_pos = projection * view * vec4(P_intersect_world, 1.0);
        float ndc_depth = clip_space_pos.z / clip_space_pos.w; // Perspective divide
        gl_FragDepth = (ndc_depth * 0.5) + 0.5; // Map from [-1, 1] to [0, 1]
    }
}