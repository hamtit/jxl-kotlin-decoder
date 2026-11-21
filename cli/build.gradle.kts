plugins {
    kotlin("jvm")
    application
}

kotlin {
    compilerOptions {
        jvmTarget.set(org.jetbrains.kotlin.gradle.dsl.JvmTarget.JVM_17)
        allWarningsAsErrors.set(providers.gradleProperty("werror").isPresent)
    }
}

java {
    sourceCompatibility = JavaVersion.VERSION_17
    targetCompatibility = JavaVersion.VERSION_17
}

dependencies {
    implementation(project(":jxl-decoder"))
}

application {
    mainClass.set("io.github.hamtit.jxl.cli.MainKt")
    applicationDefaultJvmArgs = listOf("-Xmx2g")
}
