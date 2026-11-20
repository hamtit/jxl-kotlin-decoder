plugins {
    kotlin("jvm")
    `java-library`
    `maven-publish`
}

group = "io.github.hamtit"
version = "0.1.0"

kotlin {
    // Android 12 (API 31) 以降の ART でも動くように Java 17 のバイトコードにする
    compilerOptions {
        jvmTarget.set(org.jetbrains.kotlin.gradle.dsl.JvmTarget.JVM_17)
    }
    explicitApi()
}

java {
    sourceCompatibility = JavaVersion.VERSION_17
    targetCompatibility = JavaVersion.VERSION_17
    withSourcesJar()
}

dependencies {
    testImplementation(kotlin("test"))
}

tasks.test {
    useJUnitPlatform()
    maxHeapSize = "2g"
    systemProperty("jxl.testdata", rootProject.file("testdata").absolutePath)
}

publishing {
    publications {
        create<MavenPublication>("maven") {
            artifactId = "jxl-decoder"
            from(components["java"])
        }
    }
}

